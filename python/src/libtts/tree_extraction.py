"""
Functions for extracting individual tree point clouds from a large dataset. 

The process is designed to be run in parallel for multiple trees to improve throughput on large datasets.
"""

import multiprocessing as mp
import subprocess
import os
import sys
import time
import shutil

# --- Dependency Checks ---
# Encapsulate imports in functions or check them to provide clear error messages.
try:
    # Used for reading/writing .ply files if support is added in the future.
    from plyfile import PlyData, PlyElement
except ImportError:
    print("Warning: The 'plyfile' library is not installed. Run 'pip install plyfile' to enable .ply support.")

try:
    # Used for reading .las files if a pure Python alternative to LAStools is needed.
    import laspy
except ImportError:
    pass # Not critical as LAStools is the primary method.


# --- Optional C++ Module for Object-Based Method ---
CPP_MODULE_AVAILABLE = True
try:
    from ._libtts import (tls_extract_single_trees_cpp as _tts_tls_segment,
                          generate_alpha_shape_cpp as _generate_alpha_shape
                        )
except ImportError:
    print("Info: C++ module 'libtts' not found.")
    CPP_MODULE_AVAILABLE = False


def process_single_tree(tree_id, tree_loc, clip_radius, th_alpha_sq,
                         entire_pts_file, entire_loc_file,
                         output_folder, lastools_bin_folder,
                         use_existing=False, save_intermediate=False):
    """Processes a single tree by clipping points around the detected location and running segmentation.

    Args:
        tree_id (int or str): A unique identifier for the tree.
        tree_loc (tuple[float, float]): The (x, y) coordinate of the tree's center.
        clip_radius (float): The radius in meters for the circular clip around
            the tree_loc. Suggested value: 5m.
        th_alpha_sq (float): The squared alpha threshold for the alpha shape
            generation algorithm.
        entire_pts_file (str): The file path to the complete plot-level point
            cloud (currently expects .las format).
        entire_loc_file (str): The file path to the location data used by the
            segmentation algorithm. This is typically like x,y,z,treeid, and includes
            all detected tree locations in the entire_pts_file.
        output_folder (str): The directory where the final segmented tree
            file and intermediate files will be saved.
        lastools_bin_folder (str): The file path to the 'bin' directory of your
            LAStools installation.
        use_existing (bool, optional): If True, skips processing if the final
            output file already exists. Defaults to False.
        save_intermediate (bool, optional): If True, intermediate files like the
            clipped AOI and the alpha shape are kept. Defaults to False.

    Returns:
        str | None: The path to the final segmented file if successful,
        otherwise None.

    Raises:
        RuntimeError: If the required C++ module or LAStools are not available.
        FileNotFoundError: If the input point cloud file does not exist.
        subprocess.CalledProcessError: If the `las2las` clipping command fails.
    """
    if not CPP_MODULE_AVAILABLE:
        raise RuntimeError("C++ segmentation module is not available. Cannot process tree.")

    # Check for LAStools dependency
    las2las_path = os.path.join(lastools_bin_folder, 'las2las')
    if not shutil.which(las2las_path):
        raise RuntimeError(f"LAStools 'las2las' not found or not executable at: {las2las_path}")

    if not os.path.exists(entire_pts_file):
        raise FileNotFoundError(f"Input point cloud file not found: {entire_pts_file}")

    # Define file paths using the robust os.path.join
    clipped_las_file = os.path.join(output_folder, f"tree_{tree_id}_clipped.las")
    final_seg_file = os.path.join(output_folder, f"tree_{tree_id}_a{th_alpha_sq:.3f}.pts")

    if use_existing and os.path.exists(final_seg_file):
        print(f"Skipping Tree ID {tree_id}: Final file already exists.")
        return final_seg_file

    # 1. Clip the tree's Area of Interest (AOI) from the main cloud
    try:
        if not os.path.exists(clipped_las_file):
            cmd = (
                f"{las2las_path} -i {entire_pts_file} -o {clipped_las_file} "
                f"-keep_circle {tree_loc[0]} {tree_loc[1]} {clip_radius}"
            )
            print(f"Executing for Tree ID {tree_id}: {cmd}")
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"[Error] LAStools clipping failed for Tree ID {tree_id}.")
        print(f"Stderr: {e.stderr}")
        return None

    if not os.path.exists(clipped_las_file):
        print(f"[Error] Clipped file was not created for Tree ID {tree_id}.")
        return None

    # 2. Generate alpha shape and segment the tree using the C++ module
    print(f"Running C++ segmentation for Tree ID {tree_id}...")
    try:
        as_file = _generate_alpha_shape.generate_alpha_shape(clipped_las_file, th_alpha_sq)
        seg_file_path = _tts_tls_segment.tls_extract_single_trees(
            as_file, entire_loc_file,
            th_p2trunk_distance=0.2,
            th_search_radius=0.25
        )
        # Rename the output to our standard format
        shutil.move(seg_file_path, final_seg_file)

    except Exception as e:
        print(f"[Error] C++ segmentation failed for Tree ID {tree_id}: {e}")
        return None

    # 3. Clean up intermediate files if requested
    if not save_intermediate:
        if os.path.exists(clipped_las_file):
            os.remove(clipped_las_file)
        if 'as_file' in locals() and os.path.exists(as_file):
            os.remove(as_file)

    print(f"Successfully processed Tree ID {tree_id}. Output: {final_seg_file}")
    return final_seg_file


def extract_trees_parallel(selected_tree_locs, entire_pts_file, entire_loc_file,
                           clip_radius, th_alpha_sq, output_folder,
                           lastools_bin_folder, parallel_workers=2):
    """Extracts multiple trees from a point cloud in parallel.

    Args:
        selected_tree_locs (dict): A dictionary mapping tree IDs to their
            (x, y) coordinates. Example: `{1: (x1, y1), 2: (x2, y2)}`.
        entire_pts_file (str): Path to the complete plot-level point cloud.
        entire_loc_file (str): Path to the location data file.
        clip_radius (float): The radius for the circular clip.
        th_alpha_sq (float): The squared alpha threshold for segmentation.
        output_folder (str): Directory where all outputs will be saved.
        lastools_bin_folder (str): Path to the LAStools 'bin' directory.
        parallel_workers (int, optional): The number of parallel processes to
            use. Defaults to 2.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output directory: {output_folder}")

    # Prepare a list of argument tuples for each call to process_single_tree
    tasks = []
    for tree_id, tree_loc in selected_tree_locs.items():
        tasks.append((
            tree_id, tree_loc, clip_radius, th_alpha_sq,
            entire_pts_file, entire_loc_file,
            output_folder, lastools_bin_folder,
            False, False  # Corresponds to use_existing and save_intermediate
        ))

    print(f"Created {len(tasks)} extraction tasks to process.")
    if not tasks:
        print("No tasks to run.")
        return

    # Execute tasks in parallel using a process pool
    print(f"Starting processing with {parallel_workers} parallel workers...\n")
    start_time = time.time()
    with mp.Pool(processes=parallel_workers) as pool:
        results = pool.starmap(process_single_tree, tasks)
    end_time = time.time()

    successful_tasks = [res for res in results if res is not None]
    print("\n--- Processing Complete ---")
    print(f"Successfully processed {len(successful_tasks)} out of {len(tasks)} trees.")
    print(f"Total execution time: {end_time - start_time:.2f} seconds.")
    print(f"Results saved in: {output_folder}")