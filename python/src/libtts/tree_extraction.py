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
        points (list or np.ndarray): A list of (x, y, z) tuples or numpy array of points.
        output_file (str): The path where the .ply file will be saved.
    """
    # Optimized: Use vectorized NumPy operations instead of list comprehension
    if isinstance(points, np.ndarray):
        # Already a numpy array
        if points.ndim == 2 and points.shape[1] == 3:
            # Shape: (N, 3) - already in correct format
            vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
            vertex_array = np.empty(len(points), dtype=vertex_dtype)
            vertex_array['x'] = points[:, 0].astype('f4')
            vertex_array['y'] = points[:, 1].astype('f4')
            vertex_array['z'] = points[:, 2].astype('f4')
        else:
            # Convert to list and process
            points = list(points)
            vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
            vertex_array = np.empty(len(points), dtype=vertex_dtype)
            vertex_array['x'] = np.array([p[0] for p in points], dtype='f4')
            vertex_array['y'] = np.array([p[1] for p in points], dtype='f4')
            vertex_array['z'] = np.array([p[2] for p in points], dtype='f4')
    else:
        # List of tuples - use vectorized conversion
        vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
        vertex_array = np.empty(len(points), dtype=vertex_dtype)
        # Convert to numpy array first, then assign (faster than list comprehension)
        points_arr = np.array(points, dtype='f4')
        vertex_array['x'] = points_arr[:, 0]
        vertex_array['y'] = points_arr[:, 1]
        vertex_array['z'] = points_arr[:, 2]
    
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
        np.ndarray: A NumPy array of shape (N, 3) with (x, y, z) coordinates 
        belonging to the specified tree ID. Returns empty array if no points found.
    """
    
    # Read the .ply file
    plydata = PlyData.read(seg_file)
    vertex_data = plydata['vertex'].data
    
    # Optimized: Use NumPy vectorized operations and return NumPy array directly
    # PlyData.read() returns a structured numpy array with named fields
    # Convert tree_id to int for comparison
    tree_id = int(tree_id) if isinstance(tree_id, str) else tree_id
    
    try:
        # Try structured array access (most common case)
        if hasattr(vertex_data, 'dtype') and vertex_data.dtype.names:
            # Structured array with named fields (x, y, z, label)
            print(f"using structured array access")
            labels = vertex_data['label']
            mask = labels == tree_id
            
            if np.any(mask):
                # Extract x, y, z for matching points using vectorized indexing
                x = vertex_data['x'][mask]
                y = vertex_data['y'][mask]
                z = vertex_data['z'][mask]
                # Stack into (N, 3) array - much faster than list(zip())
                points = np.column_stack([x, y, z]).astype('f4')
            else:
                points = np.empty((0, 3), dtype='f4')
        else:
            # Fallback: convert to numpy array and use positional indexing
            print(f"using positional indexing")
            vertex_arr = np.asarray(vertex_data)
            if vertex_arr.ndim == 2 and vertex_arr.shape[1] >= 4:
                # Multi-row array - columns are x, y, z, label
                labels = vertex_arr[:, 3]
                mask = labels == tree_id
                if np.any(mask):
                    # Directly extract x, y, z columns and stack - very fast
                    points = vertex_arr[mask, :3].astype('f4')
                else:
                    points = np.empty((0, 3), dtype='f4')
            else:
                # Last resort: use original method (slow, but handles edge cases)
                print(f"using original method")
                point_list = [(v[0], v[1], v[2]) for v in vertex_data if len(v) > 3 and v[3] == tree_id]
                if point_list:
                    points = np.array(point_list, dtype='f4')
                else:
                    points = np.empty((0, 3), dtype='f4')
    except (KeyError, IndexError, AttributeError):
        # If anything fails, fall back to original method (slow, but safe)
        point_list = [(v[0], v[1], v[2]) for v in vertex_data if len(v) > 3 and v[3] == tree_id]
        if point_list:
            points = np.array(point_list, dtype='f4')
        else:
            points = np.empty((0, 3), dtype='f4')
    
    return points


def process_single_tree(tree_id, tree_loc, clip_radius, th_alpha_sq,
                         entire_pts_file, entire_loc_file,
                         output_folder, lastools_bin_folder,
                         keep_random_fraction = None,
                         use_existing=False, save_intermediate=False,
                         output_target_tree=True,
                         label_propagation_method='iterative_distance_based',
                         label_propagation_kwargs=None):
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
        output_target_tree (bool, optional): If True, saves the target tree points
            to a separate file. Defaults to True.
        label_propagation_method (str, optional): The label propagation method to use. One of
            'distance_based', 'region_growing', 'region_growing_layered', 'layered_nn', or 'iterative_distance_based'.
            Defaults to 'iterative_distance_based' (recommended for dense point clouds with good coverage).
            Note: 'distance_based' may have issues with label propagation - see performance analysis.
        max_distance (float, optional): Maximum distance for distance-based method.
            Defaults to 0.1.
        search_radius (float, optional): Search radius for region growing methods.
            Defaults to 0.1.
        layer_height (float, optional): Layer height for layered methods.
            Defaults to 0.1.
        max_search_radius (float, optional): Maximum search radius for layered_nn method.
            Defaults to 0.1.
        n_jobs (int, optional): Number of parallel jobs for label propagation.
            If None, processes sequentially. Defaults to None.

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
    
    # Set default kwargs if not provided
    if label_propagation_kwargs is None:
        label_propagation_kwargs = {}

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
        # Prepare parameters for label propagation
        # Start with the output file, then merge in user-provided kwargs
        label_prop_kwargs = {}
        label_prop_kwargs.update(label_propagation_kwargs)
        
        # Generate default output filename if not provided
        if 'out_file' not in label_prop_kwargs:
            complete_lbl_file = os.path.join(output_folder, f"tree_{tree_id}_clipped_ds{keep_random_fraction:.2f}_a{th_alpha_sq:.3f}_lbl_comp.ply")
            label_prop_kwargs['out_file'] = complete_lbl_file
        else:
            complete_lbl_file = label_prop_kwargs['out_file']
        
        run_label_propagation(infile=og_clipped_las_file, labeled_file=seg_file_path, 
                              method=label_propagation_method, **label_prop_kwargs)
        ds_seg_file = seg_file_path
        seg_file_path = complete_lbl_file

    # Optionally save the target tree points
    if output_target_tree:
        import time
        save_start_time = time.time()
        target_points = get_target_tree(seg_file_path, tree_id)
        get_tree_time = time.time() - save_start_time
        
        target_file = os.path.join(output_folder, f"segtree_{tree_id}.ply")
        write_start_time = time.time()
        save_trees_as_ply(target_points, target_file)
        write_time = time.time() - write_start_time
        
        total_save_time = time.time() - save_start_time
        print(f"Saved target tree points to {target_file}")
        print(f"  Time: get_target_tree={get_tree_time:.2f}s, save_trees_as_ply={write_time:.2f}s, total={total_save_time:.2f}s")

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
                           output_target_tree=True,
                           label_propagation_method='distance_based',
                           label_propagation_kwargs=None):
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
        use_existing (bool, optional): If True, skips clipping process if
            the file already exists. Defaults to False.
        save_intermediate (bool, optional): If True, intermediate files like the
            clipped AOI and the alpha shape are kept. Defaults to False.
        output_target_tree (bool, optional): If True, saves the target tree points
            to a separate file. Defaults to True.
        label_propagation_method (str, optional): The label propagation method to use. One of
            'distance_based', 'region_growing', 'region_growing_layered', 'layered_nn', or 'iterative_distance_based'.
            Defaults to 'iterative_distance_based' (recommended for dense point clouds with good coverage).
            Note: 'distance_based' may have issues with label propagation - see performance analysis.
        label_propagation_kwargs (dict, optional): Dictionary of method-specific parameters.
            For 'distance_based': {'max_distance': 0.1, 'n_jobs': None}
            For 'region_growing_layered': {'search_radius': 0.1, 'layer_height': 0.1, 'n_jobs': None}
            For 'region_growing': {'search_radius': 0.1}
            For 'layered_nn': {'layer_height': 0.1, 'max_search_radius': 0.1}
            For 'iterative_distance_based': {'wave_distance': 0.05, 'max_iterations': 5, 
                'min_new_points': 10, 'multiple_label_strategy': 'hybrid', ...}
            Defaults to None (uses method defaults).
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output directory: {output_folder}")

    # Set default kwargs if not provided
    if label_propagation_kwargs is None:
        label_propagation_kwargs = {}
    
    # Prepare a list of argument tuples for each call to process_single_tree
    tasks = []
    for tree_id, tree_loc in selected_tree_locs.items():
        tasks.append((
            tree_id, tree_loc, clip_radius, th_alpha_sq,
            entire_pts_file, entire_loc_file,
            output_folder, lastools_bin_folder,
            keep_random_fraction,
            use_existing, save_intermediate,
            output_target_tree,
            label_propagation_method,
            label_propagation_kwargs
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