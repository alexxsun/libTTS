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

import numpy as np

from .points_downsampling import downsample_by_lastools
from .label_propagation import run_label_propagation

# --- Dependency Checks ---
# Encapsulate imports in functions or check them to provide clear error messages.
try:
    # Used for reading/writing .ply files if support is added in the future.
    from plyfile import PlyData, PlyElement
except ImportError:
    print("Warning: The 'plyfile' library is not installed. Run 'pip install plyfile' to enable .ply support.")

try:
    import laspy
except ImportError:
    print("Warning: The 'laspy' library is not installed. Run 'pip install laspy' to enable .las support.") 

# --- Optional C++ Module for Object-Based Method ---
CPP_MODULE_AVAILABLE = True
try:
    from ._libtts import (tls_extract_single_trees_cpp as _tts_tls_segment,
                          generate_alpha_shape_cpp as _generate_alpha_shape
                        )
except ImportError:
    print("Info: C++ module 'libtts' not found.")
    CPP_MODULE_AVAILABLE = False


def save_trees_as_ply(points, output_file):
    """Saves a list of points to a .ply file.

    Args:
        points (list): A list of (x, y, z) tuples representing point coordinates.
        output_file (str): The path where the .ply file will be saved.
    """
    vertex = [(p[0], p[1], p[2]) for p in points]
    vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
    vertex_array = np.array(vertex, dtype=vertex_dtype)
    ply_element = PlyElement.describe(vertex_array, 'vertex')
    PlyData([ply_element], text=False).write(output_file)
    #print(f"Saved {len(points)} points to {output_file}")
    return

def get_target_tree(seg_file, tree_id):
    """Extracts points belonging to a specific tree ID from a segmented file.
    Args:
        seg_file (str): Path to the segmented point cloud file. 
        Typically expects a .ply file, with x,y,z,label columns.
        tree_id (int or str): The tree ID to extract.
    Returns:
        list: A list of points belonging to the specified tree ID.
    """
    
    # read the .ply file
    plydata = PlyData.read(seg_file)
    vertex_data = plydata['vertex'].data
    points = [(v[0], v[1], v[2]) for v in vertex_data if v[3] == tree_id]
    return points


def process_single_tree(tree_id, tree_loc, clip_radius, th_alpha_sq,
                         entire_pts_file, entire_loc_file,
                         output_folder, lastools_bin_folder,
                         keep_random_fraction = None,
                         use_existing=False, save_intermediate=False,
                         output_target_tree=True):
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
        use_existing (bool, optional): If True, skips clipping process if
            the file already exists. Defaults to False.
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
    

    # todo: use_existing when enable, check the intermediate files for each step
    #   skip steps if they exist

    # Define file paths using the robust os.path.join
    clipped_las_file = os.path.join(output_folder, f"tree_{tree_id}_clipped.las")

    # 1. Clip the tree's Area of Interest (AOI) from the main cloud
    try:
        if not use_existing or not os.path.exists(clipped_las_file):
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

    # Optionally downsample the clipped file
    if keep_random_fraction is not None:
            # Downsample the clipped file if requested
            downsampled_file = downsample_by_lastools(infile = clipped_las_file, 
                                 lastools_bin_dir = lastools_bin_folder, keep_random_fraction = keep_random_fraction,
                                 )
            og_clipped_las_file = clipped_las_file
            clipped_las_file = downsampled_file

    # 2. Generate alpha shape and segment the tree using the C++ module
    print(f"Running C++ segmentation for Tree ID {tree_id}...")
    try:
        # Generate alpha shape
        as_file = _generate_alpha_shape(clipped_las_file, th_alpha_sq, ".ply") # ensure .ply extension
        seg_file_path = _tts_tls_segment(
            as_file, entire_loc_file,
            th_p2trunk_distance=0.2,
            th_search_radius=0.25
        )
    except Exception as e:
        print(f"[Error] C++ segmentation failed for Tree ID {tree_id}: {e}")
        return None
    
    #  Label propagation for downsampled points
    if keep_random_fraction is not None:
        complete_lbl_file = os.path.join(output_folder, f"tree_{tree_id}_clipped_ds{keep_random_fraction:.2f}_a{th_alpha_sq:.3f}_lbl_comp.ply")
        run_label_propagation(infile = og_clipped_las_file, labeled_file = seg_file_path, 
                              method='region_growing', search_radius = 0.1, out_file=complete_lbl_file)
        ds_seg_file = seg_file_path
        seg_file_path = complete_lbl_file

    # Optionally save the target tree points
    if output_target_tree:
        target_points = get_target_tree(seg_file_path, tree_id)
        target_file = os.path.join(output_folder, f"segtree_{tree_id}.ply")
        save_trees_as_ply(target_points, target_file)
        print(f"Saved target tree points to {target_file}")

    # 3. Clean up intermediate files if requested
    if not save_intermediate:

        if os.path.exists(clipped_las_file):
            os.remove(clipped_las_file)
        if 'as_file' in locals() and os.path.exists(as_file):
            os.remove(as_file)

        # todo: consider removing seg_file_path,
        # if os.path.exists(seg_file_path):
        #     os.remove(seg_file_path)
        
        
        if keep_random_fraction is not None:
            if os.path.exists(og_clipped_las_file):
                os.remove(og_clipped_las_file)
            if os.path.exists(ds_seg_file):
                os.remove(ds_seg_file)

    print(f"Successfully processed Tree ID {tree_id}. Output: {seg_file_path}")
    return seg_file_path


def extract_trees_parallel(selected_tree_locs, entire_pts_file, entire_loc_file,
                           clip_radius, th_alpha_sq, output_folder,
                           lastools_bin_folder, 
                           keep_random_fraction = None,
                           parallel_workers=2,
                           use_existing=False, save_intermediate=False,
                           output_target_tree=True):
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
            keep_random_fraction,
            use_existing, save_intermediate,
            output_target_tree 
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