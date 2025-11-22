"""Propagates labels from a small set of seed points to a larger point cloud.

This module provides two primary methods for label propagation:

1.  **Region Growing**: A method that starts from seed points and "grows"
    regions outwards, labeling neighbors within a specified radius. It processes
    points in ascending Z-height order to ensure a bottom-up propagation,
    which is useful for environments like forests.
2.  **Layered Nearest Neighbor**: A faster, more constrained method that divides
    the point cloud into horizontal layers. It labels points in each layer
    based on the nearest labeled seed within that same layer, but only if the
    seed is within a maximum search radius.

Example:
    .. code-block:: python

        import libtts

        # Using the Z-ordered region growing method
        labels = libtts.run_label_propagation(
            infile="unlabeled_cloud.ply",
            labeled_file="tree_locations.pts",
            method='region_growing',
            search_radius=0.5,
            out_file="cloud_labeled_region_growing.ply"
        )

        # Using the faster, layered nearest neighbor method
        labels_nn = libtts.run_label_propagation(
            infile="unlabeled_cloud.ply",
            labeled_file="tree_locations.pts",
            method='layered_nn',
            layer_height=1.0,
            max_search_radius=1.0,
            out_file="cloud_labeled_layered_nn.ply"
        )

"""
from typing import Optional
import numpy as np
from scipy.spatial import KDTree, cKDTree
import heapq
from plyfile import PlyData, PlyElement
import multiprocessing as mp
from multiprocessing import shared_memory
from functools import partial
import os
import time
from concurrent.futures import ThreadPoolExecutor

try:
    import laspy
except ImportError:
    print("Warning: The 'laspy' library is not installed. Run 'pip install laspy' to enable .las support.") 


# Global variables for worker processes (initialized once per worker)
_worker_kdtree = None
_worker_points = None

def _is_in_multiprocessing_context():
    """Check if we're already running inside a multiprocessing worker.
    
    Returns:
        bool: True if we're in a multiprocessing context (daemonic process)
    """
    try:
        # Check if current process is a daemonic process
        current_process = mp.current_process()
        if current_process.name != 'MainProcess':
            # We're in a worker process
            return True
        # Check if we can create a test pool (if we can't, we're likely in a worker)
        return False
    except Exception:
        # If we can't check, assume we're safe (better to try and fail gracefully)
        return False

def _init_worker(points_data):
    """Initialize worker process with shared KDTree.
    
    Args:
        points_data: The points array to build KDTree from
    """
    global _worker_kdtree, _worker_points
    _worker_points = points_data
    _worker_kdtree = KDTree(points_data)

def _process_point_parallel(args):
    """Worker function for parallel processing of points in a layer (multiprocessing).
    
    Args:
        args: Tuple of (point_idx, label, search_radius)
    
    Returns:
        Tuple of (point_idx, label, list of neighbor_indices)
    """
    global _worker_kdtree, _worker_points
    point_idx, label, search_radius = args
    
    # Use the pre-built KDTree from worker initialization
    if _worker_kdtree is None:
        raise RuntimeError("Worker not initialized. KDTree not available.")
    
    # Find all neighbors within the search radius
    neighbor_indices = _worker_kdtree.query_ball_point(_worker_points[point_idx], r=search_radius)
    
    return (point_idx, label, neighbor_indices)


def _init_worker_thread(points_data):
    """Initialize thread worker with shared KDTree.
    
    For threading, we can use the same global variables since threads share memory.
    
    Args:
        points_data: The points array to build KDTree from
    """
    global _worker_kdtree, _worker_points
    _worker_points = points_data
    _worker_kdtree = KDTree(points_data)


def _process_point_parallel_thread(args):
    """Worker function for parallel processing of points in a layer (threading).
    
    This function is used when we're in a nested multiprocessing context.
    Threads share memory, so we can use the global KDTree directly.
    
    Args:
        args: Tuple of (point_idx, label, search_radius)
    
    Returns:
        Tuple of (point_idx, label, list of neighbor_indices)
    """
    global _worker_kdtree, _worker_points
    point_idx, label, search_radius = args
    
    # Use the pre-built KDTree (shared across threads)
    if _worker_kdtree is None:
        raise RuntimeError("Thread worker not initialized. KDTree not available.")
    
    # Find all neighbors within the search radius
    # Note: KDTree.query_ball_point releases GIL, so threading works well here
    neighbor_indices = _worker_kdtree.query_ball_point(_worker_points[point_idx], r=search_radius)
    
    return (point_idx, label, neighbor_indices)


def label_points_layered_nn(
    points: np.ndarray,
    labeled_seed_points: np.ndarray,
    layer_height: float = 1.0,
    max_search_radius: float = 1.0,
    out_file: Optional[str] = None
) -> np.ndarray:
    """Labels unlabeled points using a layered Nearest Neighbor approach.

    This method processes the point cloud in horizontal layers from bottom to
    top. For each layer, it finds the nearest labeled seed point for every
    unlabeled point within that same layer. A label is only propagated if the
    seed point is within the `max_search_radius`.

    Args:
        points (np.ndarray): The target Nx3 (X,Y,Z) point cloud to be labeled.
        labeled_seed_points (np.ndarray): An Mx4 (X,Y,Z,label) array of
            initially labeled points.
        layer_height (float): The thickness of each horizontal layer for processing.
        max_search_radius (float): The maximum distance to a labeled point for
            propagation to occur.
        out_file (str, optional): If provided, saves the final labeled points
            (X,Y,Z,label) to this path. Supports .ply and .pts. Defaults to None.

    Returns:
        np.ndarray: 
            An array of shape (N,) containing the propagated labels for
            the target 'points' cloud. Unlabeled points will have a value of -1.
    """
    print(f"Starting Layered Nearest Neighbor with layer height: {layer_height} and max radius: {max_search_radius}")
    
    # Filter out seed points with invalid labels (e.g., < 1) at the beginning.
    # The label is in the 4th column (index 3).
    valid_seeds_mask = labeled_seed_points[:, 3] >= 1
    labeled_seed_points = labeled_seed_points[valid_seeds_mask]
    
    # Step 1: Create a new label array for the target points, initialized to "unlabeled" (-1).
    propagated_labels = np.full(points.shape[0], -1, dtype=int)
    
    # Get Z coordinates for both point clouds to define layers
    points_z = points[:, 2]
    seeds_z = labeled_seed_points[:, 2]
    
    # Determine the full vertical range
    min_z = min(np.min(points_z), np.min(seeds_z)) if seeds_z.size > 0 else np.min(points_z)
    max_z = max(np.max(points_z), np.max(seeds_z)) if seeds_z.size > 0 else np.max(points_z)
    
    # Step 2: Iterate through the cloud in layers from bottom to top
    for z in np.arange(min_z, max_z, layer_height):
        layer_min_z = z
        layer_max_z = z + layer_height
        
        # Find indices of points and seeds within the current layer
        layer_point_indices = np.where((points_z >= layer_min_z) & (points_z < layer_max_z))[0]
        layer_seed_indices = np.where((seeds_z >= layer_min_z) & (seeds_z < layer_max_z))[0]
        
        # If there are no seeds in this layer to propagate from, or no points to propagate to, skip
        if len(layer_seed_indices) == 0 or len(layer_point_indices) == 0:
            continue
        
        # Get the actual data for the current layer
        layer_points_to_label = points[layer_point_indices]
        layer_seeds = labeled_seed_points[layer_seed_indices]
        
        layer_seed_coords = layer_seeds[:, :3]
        layer_seed_labels = layer_seeds[:, 3].astype(int)
        
        # Build KDTree on the labeled seed points within this layer
        kdtree_layer = KDTree(layer_seed_coords)
        
        # Find nearest neighbors for the unlabeled points in this layer
        distances, indices = kdtree_layer.query(layer_points_to_label, k=1)
        
        # Create a mask for points that are within the search radius
        within_radius_mask = (distances <= max_search_radius)
        
        # If no points are close enough, skip to the next layer
        if not np.any(within_radius_mask):
            continue
        
        # Get the labels from the seeds that are close enough
        nearest_labels = layer_seed_labels[indices[within_radius_mask]]
        
        # Get the original, global indices of the points to update
        global_indices_to_update = layer_point_indices[within_radius_mask]
        
        # Assign the new labels in the main labels array
        propagated_labels[global_indices_to_update] = nearest_labels
            
    print("Layered Nearest Neighbor labeling complete.")
    # If an output file is specified, save the labels: xyzl.pts
    if out_file is not None:
        if out_file.endswith('.ply'):
            vertices = []
            for p,l in zip(points, propagated_labels):
                vertices.append((p[0], p[1], p[2], l))
            vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('label', 'i4')]
            vertex_array = np.array(vertices, dtype=vertex_dtype)
            ply_element = PlyElement.describe(vertex_array, 'vertex')
            PlyData([ply_element], text=False).write(out_file)
            print(f"Labels saved to {out_file}")
        elif out_file.endswith('.pts'):
            np.savetxt(out_file, np.column_stack((points, propagated_labels)), fmt='%.3f', delimiter=' ')
            print(f"Labels saved to {out_file}")
        else:
            print(f"Unsupported output file format: {out_file}. Labels not saved.")
    return propagated_labels


def label_points_region_growing(
    points: np.ndarray,
    labeled_seed_points: np.ndarray,
    search_radius: float = 0.5,
    out_file: Optional[str] = None
) -> np.ndarray:
    """Labels unlabeled points using a region growing method ordered by Z-height.

    This method works by starting with initial seed points and iteratively
    "growing" their labels outwards to neighboring points. A priority queue
    ensures that points with lower Z-coordinates are processed first, creating
    a bottom-up labeling effect.

    Args:
        points (np.ndarray): The target Nx3 (X,Y,Z) point cloud to be labeled.
        labeled_seed_points (np.ndarray): An Mx4 (X,Y,Z,label) array of
            initially labeled points.
        search_radius (float): The maximum distance to consider points as
            "connected" neighbors for region growing.
        out_file (str, optional): If provided, saves the final labeled points
            (X,Y,Z,label) to this path. Supports .ply and .pts. Defaults to None.

    Returns:
        np.ndarray: 
            An array of shape (N,) containing the propagated labels for
            the target 'points' cloud. Unlabeled points will have a value of -1.
    """
    print(f"Starting Z-Ordered Region Growing with radius: {search_radius}")

    # --- Step 1: Initialize Labels ---
    # Create a label array for the target points, initialized to "unlabeled" (-1).
    propagated_labels = np.full(points.shape[0], -1, dtype=int)
    
    # Filter seed points to ensure labels are valid (>= 1).
    valid_seeds_mask = labeled_seed_points[:, 3] >= 1
    valid_seeds = labeled_seed_points[valid_seeds_mask]

    if valid_seeds.shape[0] == 0:
        print("Warning: No valid seed points (label >= 1) found.")
        return propagated_labels

    # Build KDTree on the full point cloud to find where to place the seed labels.
    kdtree = KDTree(points)
    
    # Project the initial seed labels onto the nearest points in the main cloud.
    seed_coords = valid_seeds[:, :3]
    seed_labels = valid_seeds[:, 3].astype(int)
    distances, indices = kdtree.query(seed_coords, k=1)
    propagated_labels[indices] = seed_labels

    # --- Step 2: Set up Z-Ordered Priority Queue ---
    # A priority queue (min-heap) will ensure we always process lower points first.
    # The queue stores tuples of (z_coordinate, point_index).
    priority_queue = []
    
    # Add all initially labeled points to the priority queue.
    initial_labeled_indices = np.where(propagated_labels >= 1)[0]
    for idx in initial_labeled_indices:
        z_coord = points[idx, 2]
        heapq.heappush(priority_queue, (z_coord, idx))

    print(f"Found {len(np.unique(seed_labels))} unique labels to grow from.")

    # --- Step 3: Grow Regions from the Bottom Up ---
    while priority_queue:
        # Pop the point with the lowest Z-coordinate.
        current_z, current_idx = heapq.heappop(priority_queue)
        current_label = propagated_labels[current_idx]
        
        # Find all neighbors within the search radius.
        neighbor_indices = kdtree.query_ball_point(points[current_idx], r=search_radius)
        
        for neighbor_idx in neighbor_indices:
            # If the neighbor is unlabeled, label it and add it to the priority queue.
            if propagated_labels[neighbor_idx] == -1:
                propagated_labels[neighbor_idx] = current_label
                neighbor_z = points[neighbor_idx, 2]
                heapq.heappush(priority_queue, (neighbor_z, neighbor_idx))
                
    print("Z-Ordered Region Growing complete.")
        
    # If an output file is specified, save the labels: xyzl
    if out_file is not None:
        if out_file.endswith('.ply'):
            vertices = []
            for p,l in zip(points, propagated_labels):
                vertices.append((p[0], p[1], p[2], l))
            vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('label', 'i4')]
            vertex_array = np.array(vertices, dtype=vertex_dtype)
            ply_element = PlyElement.describe(vertex_array, 'vertex')
            PlyData([ply_element], text=False).write(out_file)
            print(f"Labels saved to {out_file}")
        elif out_file.endswith('.pts'):
            np.savetxt(out_file, np.column_stack((points, propagated_labels)), fmt=['%.3f', '%.3f', '%.3f', '%d'], delimiter=' ')
            print(f"Labels saved to {out_file}")
        else:
            print(f"Unsupported output file format: {out_file}. Labels not saved.")
    else:
        print("No output file specified, labels not saved.")
        
    return propagated_labels


def label_points_region_growing_layered(
    points: np.ndarray,
    labeled_seed_points: np.ndarray,
    search_radius: float = 0.5,
    layer_height: float = 0.1,
    n_jobs: Optional[int] = 4,
    out_file: Optional[str] = None
) -> np.ndarray:
    """Labels unlabeled points using a Z-layer batched region growing method.
    
    This is a parallelized version of region growing that maintains global Z-ordering
    by processing points in Z-layers sequentially, while batch-processing points
    within each layer for improved performance. This approach guarantees the same
    results as the sequential region growing algorithm while providing significant
    speedup, especially when there are many unique labels.
    
    The algorithm:
    1. Divides the point cloud into horizontal Z-layers
    2. Processes layers from bottom to top (maintains global Z-ordering)
    3. Within each layer, batch-processes all labeled points and their neighbors
    4. Uses priority queues to maintain order for newly labeled points in next layers

    Args:
        points (np.ndarray): The target Nx3 (X,Y,Z) point cloud to be labeled.
        labeled_seed_points (np.ndarray): An Mx4 (X,Y,Z,label) array of
            initially labeled points.
        search_radius (float): The maximum distance to consider points as
            "connected" neighbors for region growing.
        layer_height (float): The height of each Z-layer for batch processing.
            Smaller values maintain finer Z-ordering but create more layers.
            Defaults to 0.1.
        n_jobs (int, optional): Number of parallel workers for batch processing
            within layers. When set to a value > 1, uses multiprocessing to
            parallelize neighbor queries within each layer. Use -1 to use all
            available CPU cores. Set to None for sequential processing.
            Only used for layers with many points (>10 or 5*n_jobs). 
            Defaults to 4 (parallel processing with 4 workers).
        out_file (str, optional): If provided, saves the final labeled points
            (X,Y,Z,label) to this path. Supports .ply and .pts. Defaults to None.

    Returns:
        np.ndarray: 
            An array of shape (N,) containing the propagated labels for
            the target 'points' cloud. Unlabeled points will have a value of -1.
    """
    if n_jobs is not None and n_jobs > 1:
        print(f"Starting Z-Layer Batched Region Growing with radius: {search_radius}, layer_height: {layer_height}, n_jobs={n_jobs}")
    else:
        print(f"Starting Z-Layer Batched Region Growing with radius: {search_radius}, layer_height: {layer_height}")
    
    # --- Step 1: Initialize Labels ---
    # Create a label array for the target points, initialized to "unlabeled" (-1).
    propagated_labels = np.full(points.shape[0], -1, dtype=int)
    
    # Filter seed points to ensure labels are valid (>= 1).
    valid_seeds_mask = labeled_seed_points[:, 3] >= 1
    valid_seeds = labeled_seed_points[valid_seeds_mask]

    if valid_seeds.shape[0] == 0:
        print("Warning: No valid seed points (label >= 1) found.")
        return propagated_labels

    # Build KDTree on the full point cloud to find where to place the seed labels.
    kdtree = KDTree(points)
    
    # Project the initial seed labels onto the nearest points in the main cloud.
    seed_coords = valid_seeds[:, :3]
    seed_labels = valid_seeds[:, 3].astype(int)
    distances, indices = kdtree.query(seed_coords, k=1)
    propagated_labels[indices] = seed_labels

    # Get Z coordinates
    points_z = points[:, 2]
    
    # Determine the full vertical range
    min_z = np.min(points_z)
    max_z = np.max(points_z)
    
    # Get unique labels count
    unique_labels = len(np.unique(seed_labels))
    print(f"Found {unique_labels} unique labels to grow from.")
    print(f"Z-range: [{min_z:.3f}, {max_z:.3f}], will create ~{int((max_z - min_z) / layer_height)} layers")
    
    # --- Step 2: Set up Z-Layer Processing ---
    # We'll process layers sequentially, but batch process within each layer
    
    # Initialize priority queues for each layer
    # Use a dictionary: layer_index -> list of (z_coord, point_index, label) tuples
    layer_queues = {}
    
    # Add all initially labeled points to their respective layer queues
    initial_labeled_indices = np.where(propagated_labels >= 1)[0]
    for idx in initial_labeled_indices:
        z_coord = points_z[idx]
        layer_idx = int((z_coord - min_z) / layer_height)
        if layer_idx not in layer_queues:
            layer_queues[layer_idx] = []
        layer_queues[layer_idx].append((z_coord, idx, propagated_labels[idx]))
    
    # Sort each layer's queue by Z-coordinate (maintains ordering within layer)
    for layer_idx in layer_queues:
        layer_queues[layer_idx].sort(key=lambda x: x[0])
    
    # --- Step 3: Process Layers from Bottom to Top ---
    # Always process the layer with the lowest Z that has points (maintains global Z-ordering)
    processed_layers = 0
    total_points_labeled = len(initial_labeled_indices)
    
    # Keep processing until no layers have points
    while layer_queues:
        # Find the layer with the lowest Z that has points
        active_layers = [idx for idx in layer_queues.keys() if len(layer_queues[idx]) > 0]
        if not active_layers:
            break
        
        # Process the lowest layer (maintains global Z-ordering)
        current_layer_idx = min(active_layers)
        current_layer_queue = layer_queues[current_layer_idx]
        next_layer_points = {}  # layer_idx -> list of (z, idx, label)
        
        # Determine if we should use parallel processing
        # Check if we're already in a multiprocessing context
        in_mp_context = _is_in_multiprocessing_context()
        use_parallel = (n_jobs is not None and n_jobs > 1 and 
                       len(current_layer_queue) > max(10, n_jobs * 5))
        
        if use_parallel:
            # Parallel processing: batch query neighbors for all points in layer
            # Determine number of workers
            num_workers = n_jobs if n_jobs > 0 else mp.cpu_count()
            num_workers = min(num_workers, len(current_layer_queue), mp.cpu_count())
            
            # Prepare arguments for parallel processing
            # We process points in Z-order, so maintain the queue order
            process_args = [(idx, label, search_radius) 
                          for _, idx, label in current_layer_queue]
            
            # Hybrid approach: Use threading if in nested context, multiprocessing otherwise
            parallel_success = False
            if in_mp_context:
                # We're in a nested multiprocessing context - use threading
                #print("Using threading for parallel processing - nested context detected")
                # Initialize the global KDTree for threads (shared memory)
                _init_worker_thread(points)
                
                try:
                    # Use ThreadPoolExecutor (works in daemonic processes)
                    #print(f"  Using threading for parallel processing ({num_workers} workers) - nested context detected")
                    with ThreadPoolExecutor(max_workers=num_workers) as executor:
                        results = list(executor.map(_process_point_parallel_thread, process_args))
                    parallel_success = True
                except Exception as e:
                    # If threading fails, fall back to sequential
                    print(f"Warning: Threading failed, falling back to sequential: {e}")
                    use_parallel = False
            else:
                # We're in the main process - use multiprocessing (better performance)
                #print(f"Using multiprocessing for parallel processing - main process")
                try:
                    # Process points in parallel with initialized workers
                    #print(f"  Using multiprocessing for parallel processing ({num_workers} workers) - main process")
                    with mp.Pool(processes=num_workers, initializer=_init_worker, initargs=(points,)) as pool:
                        results = pool.map(_process_point_parallel, process_args)
                    parallel_success = True
                except (RuntimeError, AssertionError) as e:
                    # If we can't create a pool, fall back to sequential
                    if "daemonic" in str(e) or "not allowed to have children" in str(e):
                        # This shouldn't happen if detection is correct, but handle it anyway
                        print(f"Warning: Multiprocessing failed, falling back to sequential: {e}")
                        use_parallel = False
                    else:
                        raise
            
            if parallel_success:
                # Process results in Z-order (maintains ordering)
                for i, (z_coord, original_idx, original_label) in enumerate(current_layer_queue):
                    # Skip if this point was already processed by a different label
                    if propagated_labels[original_idx] != original_label:
                        continue
                    
                    # Get the parallel query result
                    result_idx, result_label, neighbor_indices = results[i]
                    
                    # Verify we got the right result
                    if result_idx != original_idx:
                        # Fallback: query directly if mismatch (shouldn't happen)
                        neighbor_indices = kdtree.query_ball_point(points[original_idx], r=search_radius)
                    
                    # Vectorized filtering of unlabeled neighbors
                    neighbor_array = np.array(neighbor_indices, dtype=int)
                    if len(neighbor_array) > 0:
                        unlabeled_mask = propagated_labels[neighbor_array] == -1
                        unlabeled_neighbors = neighbor_array[unlabeled_mask]
                        
                        # Vectorized label assignment
                        if len(unlabeled_neighbors) > 0:
                            propagated_labels[unlabeled_neighbors] = original_label
                            
                            # Add to appropriate layer queues
                            neighbor_z_coords = points_z[unlabeled_neighbors]
                            for neighbor_idx, neighbor_z in zip(unlabeled_neighbors, neighbor_z_coords):
                                neighbor_layer_idx = int((neighbor_z - min_z) / layer_height)
                                
                                if neighbor_layer_idx not in next_layer_points:
                                    next_layer_points[neighbor_layer_idx] = []
                                next_layer_points[neighbor_layer_idx].append((neighbor_z, neighbor_idx, original_label))
                                
                                total_points_labeled += 1
                
                # Clear the current layer queue since we processed all points
                current_layer_queue.clear()
        
        if not use_parallel:
            # Sequential processing (original implementation)
            while current_layer_queue:
                # Pop the point with lowest Z in this layer
                current_z, current_idx, current_label = current_layer_queue.pop(0)
                
                # Skip if this point was already processed by a different label
                # (can happen if point was added to queue multiple times)
                if propagated_labels[current_idx] != current_label:
                    continue
                
                # Find all neighbors within the search radius
                neighbor_indices = kdtree.query_ball_point(points[current_idx], r=search_radius)
                
                # Batch process neighbors: filter unlabeled ones
                unlabeled_neighbors = []
                for neighbor_idx in neighbor_indices:
                    if propagated_labels[neighbor_idx] == -1:
                        unlabeled_neighbors.append(neighbor_idx)
                
                # Label all unlabeled neighbors and add them to appropriate layer queues
                for neighbor_idx in unlabeled_neighbors:
                    propagated_labels[neighbor_idx] = current_label
                    neighbor_z = points_z[neighbor_idx]
                    
                    # Determine which layer this neighbor belongs to
                    neighbor_layer_idx = int((neighbor_z - min_z) / layer_height)
                    
                    # Add to the appropriate layer queue
                    if neighbor_layer_idx not in next_layer_points:
                        next_layer_points[neighbor_layer_idx] = []
                    next_layer_points[neighbor_layer_idx].append((neighbor_z, neighbor_idx, current_label))
                    
                    total_points_labeled += 1
        
        # Remove current layer from queues if empty
        if len(layer_queues[current_layer_idx]) == 0:
            del layer_queues[current_layer_idx]
        
        # Merge newly labeled points into layer queues
        for next_layer_idx, new_points in next_layer_points.items():
            if next_layer_idx not in layer_queues:
                layer_queues[next_layer_idx] = []
            layer_queues[next_layer_idx].extend(new_points)
            # Sort by Z-coordinate to maintain ordering
            layer_queues[next_layer_idx].sort(key=lambda x: x[0])
        
        processed_layers += 1
    
    print(f"Z-Layer Batched Region Growing complete. Processed {processed_layers} layers, labeled {total_points_labeled} points.")
    
    # If an output file is specified, save the labels: xyzl
    if out_file is not None:
        if out_file.endswith('.ply'):
            # Optimized: Create structured array directly from numpy arrays (vectorized, no Python loop)
            vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('label', 'i4')]
            vertex_array = np.empty(len(points), dtype=vertex_dtype)
            vertex_array['x'] = points[:, 0].astype('f4')
            vertex_array['y'] = points[:, 1].astype('f4')
            vertex_array['z'] = points[:, 2].astype('f4')
            vertex_array['label'] = propagated_labels.astype('i4')
            ply_element = PlyElement.describe(vertex_array, 'vertex')
            PlyData([ply_element], text=False).write(out_file)
            print(f"Labels saved to {out_file}")
        elif out_file.endswith('.pts'):
            np.savetxt(out_file, np.column_stack((points, propagated_labels)), 
                      fmt=['%.3f', '%.3f', '%.3f', '%d'], delimiter=' ')
            print(f"Labels saved to {out_file}")
        else:
            print(f"Unsupported output file format: {out_file}. Labels not saved.")
    else:
        print("No output file specified, labels not saved.")
        
    return propagated_labels


def label_points_distance_based(
    points: np.ndarray,
    labeled_seed_points: np.ndarray,
    max_distance: float = 0.25,
    n_jobs: Optional[int] = None,
    batch_size: int = 10000,
    out_file: Optional[str] = None
) -> np.ndarray:
    """Label points using distance-based nearest neighbor assignment.
    
    This method builds a KDTree on labeled seed points only (much smaller than
    the full point cloud), then for each unlabeled point, finds the nearest
    labeled point and assigns its label if within the maximum distance threshold.
    
    This method is optimized for dense point clouds where the number of labeled
    seed points is much smaller than the total number of points.
    
    Args:
        points (np.ndarray): The target Nx3 (X,Y,Z) point cloud to be labeled.
        labeled_seed_points (np.ndarray): An Mx4 (X,Y,Z,label) array of
            initially labeled points.
        max_distance (float): Maximum distance to nearest labeled point for
            assignment. Points beyond this distance remain unlabeled. Defaults to 0.25m.
        n_jobs (int, optional): Number of parallel jobs for batch processing.
            If None, processes sequentially. If > 1, uses multiprocessing.
            Defaults to None.
        batch_size (int): Number of points to process in each batch when using
            parallel processing. Defaults to 10000.
        out_file (str, optional): If provided, saves the final labeled points
            (X,Y,Z,label) to this path. Supports .ply and .pts. Defaults to None.
    
    Returns:
        np.ndarray: 
            An array of shape (N,) containing the propagated labels for
            the target 'points' cloud. Unlabeled points will have a value of -1.
    """
    print(f"Starting Distance-Based Label Propagation with max_distance: {max_distance}")
    
    # Filter out seed points with invalid labels (e.g., < 1)
    valid_seeds_mask = labeled_seed_points[:, 3] >= 1
    labeled_seed_points = labeled_seed_points[valid_seeds_mask]
    
    if labeled_seed_points.shape[0] == 0:
        print("Warning: No valid labeled seed points found. Returning all unlabeled.")
        return np.full(points.shape[0], -1, dtype=int)
    
    print(f"Using {labeled_seed_points.shape[0]} labeled seed points for KDTree")
    
    # Initialize label array
    propagated_labels = np.full(points.shape[0], -1, dtype=int)
    
    # Extract labeled point coordinates and labels
    labeled_points = labeled_seed_points[:, :3]  # X, Y, Z
    # Convert labels to integers (handle float labels like 1.0 -> 1)
    labeled_point_labels = np.round(labeled_seed_points[:, 3]).astype(int)  # labels
    
    # Step 1: Project seed labels onto the full point cloud
    # This is necessary because seed points might come from a downsampled file
    # and may not exactly match points in the full cloud
    print("Projecting seed labels onto full point cloud...")
    full_cloud_kdtree = cKDTree(points)
    seed_distances, seed_indices = full_cloud_kdtree.query(labeled_points, k=1)
    
    # Project labels within a small tolerance (e.g., 0.01m for exact matches, or use max_distance for fuzzy matching)
    projection_tolerance = min(0.01, max_distance * 0.1)  # Use 10% of max_distance or 0.01m, whichever is smaller
    projected_mask = seed_distances <= projection_tolerance
    num_projected = projected_mask.sum()
    
    if num_projected > 0:
        propagated_labels[seed_indices[projected_mask]] = labeled_point_labels[projected_mask]
        print(f"Projected {num_projected} seed labels onto full cloud (tolerance: {projection_tolerance}m)")
    else:
        # If no exact matches, use a larger tolerance (up to max_distance)
        print(f"Warning: No seed points matched within {projection_tolerance}m. Trying with max_distance tolerance...")
        projected_mask = seed_distances <= max_distance
        num_projected = projected_mask.sum()
        if num_projected > 0:
            propagated_labels[seed_indices[projected_mask]] = labeled_point_labels[projected_mask]
            print(f"Projected {num_projected} seed labels onto full cloud (tolerance: {max_distance}m)")
        else:
            print(f"Warning: No seed points matched even with {max_distance}m tolerance. "
                  f"Min distance: {seed_distances.min():.4f}m, Max distance: {seed_distances.max():.4f}m")
    
    # Step 2: Build KDTree on currently labeled points in the full cloud for propagation
    # Find points that are now labeled after projection
    labeled_mask = propagated_labels >= 1
    if not np.any(labeled_mask):
        print("Warning: No points were labeled after projection. Cannot propagate.")
        return propagated_labels
    
    labeled_points_in_cloud = points[labeled_mask]
    labeled_labels_in_cloud = propagated_labels[labeled_mask]
    
    print(f"Building KDTree on {len(labeled_points_in_cloud)} labeled points in full cloud...")
    kdtree = cKDTree(labeled_points_in_cloud)  # Use cKDTree for better performance
    
    # Find unlabeled points (points that don't already have labels)
    unlabeled_mask = propagated_labels < 1
    
    # Find indices of unlabeled points
    unlabeled_indices = np.where(unlabeled_mask)[0]
    unlabeled_points = points[unlabeled_indices]
    
    print(f"Labeling {len(unlabeled_indices)} unlabeled points...")
    
    # Process points in batches for better performance and memory management
    if n_jobs is not None and n_jobs > 1 and not _is_in_multiprocessing_context():
        # Use multiprocessing for parallel batch processing
        print(f"Using multiprocessing with {n_jobs} jobs...")
        try:
            with mp.Pool(processes=n_jobs, initializer=_init_worker_distance_based, 
                       initargs=(labeled_points_in_cloud, labeled_labels_in_cloud)) as pool:
                # Split unlabeled points into batches
                num_batches = (len(unlabeled_indices) + batch_size - 1) // batch_size
                batch_results = []
                
                for i in range(num_batches):
                    start_idx = i * batch_size
                    end_idx = min((i + 1) * batch_size, len(unlabeled_indices))
                    batch_indices = unlabeled_indices[start_idx:end_idx]
                    batch_points = unlabeled_points[start_idx:end_idx]
                    
                    batch_results.append(pool.apply_async(
                        _process_batch_distance_based,
                        (batch_points, batch_indices, max_distance)
                    ))
                
                # Collect results
                for result in batch_results:
                    batch_indices, batch_labels = result.get()
                    propagated_labels[batch_indices] = batch_labels
        except (RuntimeError, AssertionError) as e:
            print(f"Warning: Could not create multiprocessing pool ({e}). Falling back to sequential processing.")
            # Fall back to sequential processing
            n_jobs = None
    
    if n_jobs is None or n_jobs == 1:
        # Sequential processing (or fallback from failed multiprocessing)
        print("Processing points sequentially...")
        # Query all unlabeled points at once (cKDTree supports batch queries)
        # Note: workers parameter may not be available in all scipy versions
        try:
            distances, indices = kdtree.query(unlabeled_points, k=1, workers=1)
        except TypeError:
            # Fallback for older scipy versions
            distances, indices = kdtree.query(unlabeled_points, k=1)
        
        # Assign labels if within threshold
        within_threshold = distances <= max_distance
        propagated_labels[unlabeled_indices[within_threshold]] = labeled_labels_in_cloud[indices[within_threshold]]
    
    # Count labeled points
    num_labeled = np.sum(propagated_labels >= 1)
    print(f"Distance-Based Label Propagation complete. Labeled {num_labeled} points out of {points.shape[0]} total.")
    
    # If an output file is specified, save the labels
    if out_file is not None:
        write_start_time = time.time()
        if out_file.endswith('.ply'):
            # Prepare vertex data (optimized: vectorized, no Python loop)
            prep_start = time.time()
            vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('label', 'i4')]
            vertex_array = np.empty(len(points), dtype=vertex_dtype)
            vertex_array['x'] = points[:, 0].astype('f4')
            vertex_array['y'] = points[:, 1].astype('f4')
            vertex_array['z'] = points[:, 2].astype('f4')
            vertex_array['label'] = propagated_labels.astype('i4')
            prep_time = time.time() - prep_start
            
            # Write to file
            write_file_start = time.time()
            ply_element = PlyElement.describe(vertex_array, 'vertex')
            PlyData([ply_element], text=False).write(out_file)
            write_file_time = time.time() - write_file_start
            
            total_write_time = time.time() - write_start_time
            print(f"Labels saved to {out_file}")
            print(f"  Write time: prep={prep_time:.2f}s, write_file={write_file_time:.2f}s, total={total_write_time:.2f}s")
        elif out_file.endswith('.pts'):
            # Prepare data
            prep_start = time.time()
            data = np.column_stack((points, propagated_labels))
            prep_time = time.time() - prep_start
            
            # Write to file
            write_file_start = time.time()
            np.savetxt(out_file, data, fmt=['%.3f', '%.3f', '%.3f', '%d'], delimiter=' ')
            write_file_time = time.time() - write_file_start
            
            total_write_time = time.time() - write_start_time
            print(f"Labels saved to {out_file}")
            print(f"  Write time: prep={prep_time:.2f}s, write_file={write_file_time:.2f}s, total={total_write_time:.2f}s")
        else:
            print(f"Unsupported output file format: {out_file}. Labels not saved.")
    else:
        print("No output file specified, labels not saved.")
    
    return propagated_labels


# Worker initialization for distance-based parallel processing
_worker_labeled_points = None
_worker_labeled_labels = None
_worker_labeled_kdtree = None

def _init_worker_distance_based(labeled_points, labeled_labels):
    """Initialize worker process with labeled points KDTree."""
    global _worker_labeled_points, _worker_labeled_labels, _worker_labeled_kdtree
    _worker_labeled_points = labeled_points
    _worker_labeled_labels = labeled_labels
    _worker_labeled_kdtree = cKDTree(labeled_points)

def _process_batch_distance_based(batch_points, batch_indices, max_distance):
    """Process a batch of points for distance-based labeling (multiprocessing worker)."""
    global _worker_labeled_kdtree, _worker_labeled_labels
    distances, indices = _worker_labeled_kdtree.query(batch_points, k=1)
    within_threshold = distances <= max_distance
    batch_labels = np.full(len(batch_points), -1, dtype=int)
    batch_labels[within_threshold] = _worker_labeled_labels[indices[within_threshold]]
    return batch_indices, batch_labels


def _apply_hybrid_strategy(distances, indices, labeled_labels, wave_distance, weight_tolerance=1e-6):
    """Apply hybrid strategy: majority → weighted → closest.
    
    Args:
        distances: Array of distances to labeled points (shape: N × k)
        indices: Array of indices to labeled points (shape: N × k)
        labeled_labels: Array of labels for labeled points
        wave_distance: Maximum distance threshold
        weight_tolerance: Tolerance for considering weights as equal
    
    Returns:
        Array of assigned labels (shape: N,), -1 for unlabeled
    """
    n_points = len(distances)
    assigned_labels = np.full(n_points, -1, dtype=int)
    
    for i in range(n_points):
        # Find points within threshold
        within_threshold = distances[i] <= wave_distance
        if not np.any(within_threshold):
            continue  # No labeled points within threshold
        
        nearby_labels = labeled_labels[indices[i][within_threshold]]
        nearby_distances = distances[i][within_threshold]
        
        # Step 1: Try majority
        unique_labels, counts = np.unique(nearby_labels, return_counts=True)
        max_count = np.max(counts)
        max_count_indices = np.where(counts == max_count)[0]
        
        if len(max_count_indices) == 1:
            # Clear majority - use it
            assigned_labels[i] = unique_labels[max_count_indices[0]]
            continue
        
        # Step 2: No clear majority, try weighted
        label_weights = {}
        for label, dist in zip(nearby_labels, nearby_distances):
            weight = 1.0 / (dist + 1e-6)  # Inverse distance weight
            label_weights[label] = label_weights.get(label, 0.0) + weight
        
        # Find label with maximum weight
        max_weight = max(label_weights.values())
        max_weight_labels = [label for label, weight in label_weights.items() 
                            if abs(weight - max_weight) < weight_tolerance]
        
        if len(max_weight_labels) == 1:
            # Clear winner by weight - use it
            assigned_labels[i] = max_weight_labels[0]
            continue
        
        # Step 3: Weights are tied (or no clear winner), use closest
        # Find the closest point among those within threshold
        min_dist_idx = np.argmin(nearby_distances)
        assigned_labels[i] = nearby_labels[min_dist_idx]
    
    return assigned_labels


def _process_batch_majority_shared(batch_info):
    """Process a batch of points using vectorized operations with shared memory.
    
    This is a worker function for parallel processing of the majority strategy.
    It processes a batch of points using fully vectorized NumPy operations.
    
    Args:
        batch_info: Tuple containing:
            - start_idx: Start index of batch
            - end_idx: End index of batch
            - shm_distances_name: Shared memory name for distances array
            - shm_indices_name: Shared memory name for indices array
            - shm_labels_name: Shared memory name for labeled_labels array
            - distances_shape: Shape of distances array
            - indices_shape: Shape of indices array
            - labeled_labels_shape: Shape of labeled_labels array
            - wave_distance: Distance threshold
    
    Returns:
        Tuple of (start_idx, end_idx, batch_labels) where batch_labels is the
        assigned labels for this batch.
    """
    (start_idx, end_idx, shm_distances_name, shm_indices_name, shm_labels_name,
     distances_shape, indices_shape, labeled_labels_shape, wave_distance) = batch_info
    
    # Attach to shared memory
    shm_dist = shared_memory.SharedMemory(name=shm_distances_name)
    shm_idx = shared_memory.SharedMemory(name=shm_indices_name)
    shm_lbl = shared_memory.SharedMemory(name=shm_labels_name)
    
    try:
        distances = np.ndarray(distances_shape, dtype=np.float64, buffer=shm_dist.buf)
        indices = np.ndarray(indices_shape, dtype=np.int64, buffer=shm_idx.buf)
        labeled_labels = np.ndarray(labeled_labels_shape, dtype=np.int32, buffer=shm_lbl.buf)
        
        batch_size = end_idx - start_idx
        batch_labels = np.full(batch_size, -1, dtype=np.int32)
        
        # Extract batch data
        batch_distances = distances[start_idx:end_idx]  # Shape: (batch_size, k)
        batch_indices = indices[start_idx:end_idx]     # Shape: (batch_size, k)
        
        # Vectorized: create mask for neighbors within threshold
        valid_mask = batch_distances <= wave_distance  # Shape: (batch_size, k)
        
        # Vectorized: count valid neighbors per point
        valid_neighbor_counts = np.sum(valid_mask, axis=1)  # Shape: (batch_size,)
        
        # Process points with 1 neighbor (vectorized - fast path)
        single_mask = valid_neighbor_counts == 1
        if np.any(single_mask):
            # Find first valid neighbor for each point
            first_valid_idx = np.argmax(valid_mask[single_mask], axis=1)
            single_point_indices = np.where(single_mask)[0]
            neighbor_indices_single = batch_indices[single_mask, first_valid_idx]
            
            # Safety check: filter invalid indices
            max_valid_idx = len(labeled_labels) - 1
            valid_single_mask = neighbor_indices_single <= max_valid_idx
            if np.any(valid_single_mask):
                batch_labels[single_point_indices[valid_single_mask]] = labeled_labels[neighbor_indices_single[valid_single_mask]]
        
        # Process points with 2+ neighbors (vectorized grouping)
        multiple_mask = valid_neighbor_counts >= 2
        if np.any(multiple_mask):
            k = batch_distances.shape[1]
            N_batch = batch_size
            
            # Create point indices for each neighbor: [0,0,...,0, 1,1,...,1, ...]
            point_indices = np.repeat(np.arange(N_batch), k)  # Shape: (N_batch*k,)
            
            # Flatten valid mask and get valid point-neighbor pairs
            valid_flat = valid_mask.ravel()  # Shape: (N_batch*k,)
            
            if np.any(valid_flat):
                # Get corresponding neighbor indices and labels (vectorized)
                neighbor_indices_flat = batch_indices.ravel()[valid_flat]
                valid_point_flat = point_indices[valid_flat]
                
                # Filter to only points with multiple neighbors
                multiple_point_indices = np.where(multiple_mask)[0]
                multiple_neighbors_point_mask = np.isin(valid_point_flat, multiple_point_indices)
                neighbor_indices_flat = neighbor_indices_flat[multiple_neighbors_point_mask]
                valid_point_flat = valid_point_flat[multiple_neighbors_point_mask]
                
                # Safety check: filter out any invalid indices
                max_valid_idx = len(labeled_labels) - 1
                valid_idx_mask = neighbor_indices_flat <= max_valid_idx
                neighbor_indices_flat = neighbor_indices_flat[valid_idx_mask]
                valid_point_flat = valid_point_flat[valid_idx_mask]
                
                if len(neighbor_indices_flat) > 0:
                    neighbor_labels_flat = labeled_labels[neighbor_indices_flat]
                    
                    # Group by point and find majority label for each group
                    # Sort by point index for efficient grouping
                    sort_idx = np.argsort(valid_point_flat)
                    sorted_points = valid_point_flat[sort_idx]
                    sorted_labels = neighbor_labels_flat[sort_idx]
                    
                    # Find group boundaries
                    unique_points, group_starts = np.unique(sorted_points, return_index=True)
                    group_ends = np.append(group_starts[1:], len(sorted_points))
                    
                    # Process each group (vectorized bincount within each group)
                    for i, point_idx in enumerate(unique_points):
                        start = group_starts[i]
                        end = group_ends[i]
                        point_labels = sorted_labels[start:end]
                        
                        if len(point_labels) > 0:
                            # Count labels using vectorized bincount
                            max_label = point_labels.max()
                            if max_label >= 0:
                                counts = np.bincount(point_labels, minlength=max_label + 1)
                                majority_label = np.argmax(counts)
                                batch_labels[point_idx] = majority_label
        
        return start_idx, end_idx, batch_labels
    
    finally:
        # Clean up shared memory attachments (don't unlink - main process will do that)
        shm_dist.close()
        shm_idx.close()
        shm_lbl.close()


def _process_batch_majority_thread(batch_info):
    """Process a batch of points using vectorized operations (threading version).
    
    This is a worker function for thread-based parallel processing of the majority strategy.
    Since threads share memory, we can pass arrays directly instead of using shared_memory.
    
    Args:
        batch_info: Tuple containing:
            - start_idx: Start index of batch
            - end_idx: End index of batch
            - distances: Distances array (shared by threads)
            - indices: Indices array (shared by threads)
            - labeled_labels: Labeled labels array (shared by threads)
            - wave_distance: Distance threshold
    
    Returns:
        Tuple of (start_idx, end_idx, batch_labels) where batch_labels is the
        assigned labels for this batch.
    """
    (start_idx, end_idx, distances, indices, labeled_labels, wave_distance) = batch_info
    
    batch_size = end_idx - start_idx
    batch_labels = np.full(batch_size, -1, dtype=np.int32)
    
    # Extract batch data
    batch_distances = distances[start_idx:end_idx]  # Shape: (batch_size, k)
    batch_indices = indices[start_idx:end_idx]     # Shape: (batch_size, k)
    
    # Vectorized: create mask for neighbors within threshold
    valid_mask = batch_distances <= wave_distance  # Shape: (batch_size, k)
    
    # Vectorized: count valid neighbors per point
    valid_neighbor_counts = np.sum(valid_mask, axis=1)  # Shape: (batch_size,)
    
    # Process points with 1 neighbor (vectorized - fast path)
    single_mask = valid_neighbor_counts == 1
    if np.any(single_mask):
        # Find first valid neighbor for each point
        first_valid_idx = np.argmax(valid_mask[single_mask], axis=1)
        single_point_indices = np.where(single_mask)[0]
        neighbor_indices_single = batch_indices[single_mask, first_valid_idx]
        
        # Safety check: filter invalid indices
        max_valid_idx = len(labeled_labels) - 1
        valid_single_mask = neighbor_indices_single <= max_valid_idx
        if np.any(valid_single_mask):
            batch_labels[single_point_indices[valid_single_mask]] = labeled_labels[neighbor_indices_single[valid_single_mask]]
    
    # Process points with 2+ neighbors (vectorized grouping)
    multiple_mask = valid_neighbor_counts >= 2
    if np.any(multiple_mask):
        k = batch_distances.shape[1]
        N_batch = batch_size
        
        # Create point indices for each neighbor: [0,0,...,0, 1,1,...,1, ...]
        point_indices = np.repeat(np.arange(N_batch), k)  # Shape: (N_batch*k,)
        
        # Flatten valid mask and get valid point-neighbor pairs
        valid_flat = valid_mask.ravel()  # Shape: (N_batch*k,)
        
        if np.any(valid_flat):
            # Get corresponding neighbor indices and labels (vectorized)
            neighbor_indices_flat = batch_indices.ravel()[valid_flat]
            valid_point_flat = point_indices[valid_flat]
            
            # Filter to only points with multiple neighbors
            multiple_point_indices = np.where(multiple_mask)[0]
            multiple_neighbors_point_mask = np.isin(valid_point_flat, multiple_point_indices)
            neighbor_indices_flat = neighbor_indices_flat[multiple_neighbors_point_mask]
            valid_point_flat = valid_point_flat[multiple_neighbors_point_mask]
            
            # Safety check: filter out any invalid indices
            max_valid_idx = len(labeled_labels) - 1
            valid_idx_mask = neighbor_indices_flat <= max_valid_idx
            neighbor_indices_flat = neighbor_indices_flat[valid_idx_mask]
            valid_point_flat = valid_point_flat[valid_idx_mask]
            
            if len(neighbor_indices_flat) > 0:
                neighbor_labels_flat = labeled_labels[neighbor_indices_flat]
                
                # Group by point and find majority label for each group
                # Sort by point index for efficient grouping
                sort_idx = np.argsort(valid_point_flat)
                sorted_points = valid_point_flat[sort_idx]
                sorted_labels = neighbor_labels_flat[sort_idx]
                
                # Find group boundaries
                unique_points, group_starts = np.unique(sorted_points, return_index=True)
                group_ends = np.append(group_starts[1:], len(sorted_points))
                
                # Process each group (vectorized bincount within each group)
                for i, point_idx in enumerate(unique_points):
                    start = group_starts[i]
                    end = group_ends[i]
                    point_labels = sorted_labels[start:end]
                    
                    if len(point_labels) > 0:
                        # Count labels using vectorized bincount
                        max_label = point_labels.max()
                        if max_label >= 0:
                            counts = np.bincount(point_labels, minlength=max_label + 1)
                            majority_label = np.argmax(counts)
                            batch_labels[point_idx] = majority_label
    
    return start_idx, end_idx, batch_labels


def label_points_iterative_distance_based(
    points: np.ndarray,
    labeled_seed_points: np.ndarray,
    wave_distance: float = 0.05,
    max_iterations: int = 5,
    min_new_points: int = 10,
    multiple_label_strategy: str = 'closest',
    weight_tolerance: float = 1e-6,
    out_file: Optional[str] = None,
    n_jobs: Optional[int] = 8
) -> np.ndarray:
    """Label points using iterative distance-based wave propagation.
    
    This method propagates labels in waves, where each wave expands labels by a small
    distance (wave_distance) from currently labeled points. This creates gradual expansion
    similar to region growing but using distance queries.
    
    Args:
        points (np.ndarray): The target Nx3 (X,Y,Z) point cloud to be labeled.
        labeled_seed_points (np.ndarray): An Mx4 (X,Y,Z,label) array of
            initially labeled points.
        wave_distance (float): Distance to propagate in each wave. Defaults to 0.05.
        max_iterations (int): Maximum number of waves. Defaults to 5.
        min_new_points (int): Stop if fewer than this many points are labeled in a wave.
            Defaults to 10.
        multiple_label_strategy (str): Strategy for handling multiple labeled points.
            One of 'closest', 'hybrid', 'majority', 'same_only', 'weighted'.
            Defaults to 'closest' (recommended - fastest and good quality).
            Note: Other strategies ('majority', 'hybrid', 'weighted') are significantly slower
            and not recommended until further optimization.
        weight_tolerance (float): Tolerance for considering weights as equal in hybrid
            strategy. Defaults to 1e-6.
        out_file (str, optional): If provided, saves the final labeled points
            (X,Y,Z,label) to this path. Supports .ply and .pts. Defaults to None.
        n_jobs (int, optional): Number of parallel workers for 'majority' strategy.
            Uses shared memory for read-only data. Defaults to 8. Set to None or 1
            for sequential processing.
    
    Returns:
        np.ndarray: 
            An array of shape (N,) containing the propagated labels for
            the target 'points' cloud. Unlabeled points will have a value of -1.
    """
    print(f"Starting Iterative Distance-Based Wave Propagation")
    print(f"  wave_distance: {wave_distance}m, max_iterations: {max_iterations}, "
          f"min_new_points: {min_new_points}, strategy: {multiple_label_strategy}")
    
    # Filter out seed points with invalid labels
    valid_seeds_mask = labeled_seed_points[:, 3] >= 1
    labeled_seed_points = labeled_seed_points[valid_seeds_mask]
    
    if labeled_seed_points.shape[0] == 0:
        print("Warning: No valid labeled seed points found. Returning all unlabeled.")
        return np.full(points.shape[0], -1, dtype=int)
    
    print(f"Using {labeled_seed_points.shape[0]} labeled seed points")
    
    # Initialize label array
    propagated_labels = np.full(points.shape[0], -1, dtype=int)
    
    # Initialize with seed labels
    # Project seed labels onto nearest points in main cloud
    seed_points = labeled_seed_points[:, :3]
    seed_labels = np.round(labeled_seed_points[:, 3]).astype(int)
    
    # Build KDTree on main point cloud to find nearest points for seeds
    main_kdtree = cKDTree(points)
    seed_distances, seed_indices = main_kdtree.query(seed_points, k=1)
    
    # Assign seed labels to nearest points (within a reasonable tolerance)
    # Use a tolerance that accounts for potential coordinate differences from downsampling
    # Start with a small tolerance, but fall back to larger if needed
    seed_tolerance = min(0.01, wave_distance)  # Use 0.01m or wave_distance, whichever is smaller
    
    projected_mask = seed_distances <= seed_tolerance
    num_projected = projected_mask.sum()
    
    if num_projected > 0:
        propagated_labels[seed_indices[projected_mask]] = seed_labels[projected_mask]
        print(f"Projected {num_projected} seed labels onto full cloud (tolerance: {seed_tolerance}m)")
    else:
        # If no matches with small tolerance, try with wave_distance
        print(f"Warning: No seed points matched within {seed_tolerance}m. Trying with wave_distance tolerance...")
        seed_tolerance = wave_distance
        projected_mask = seed_distances <= seed_tolerance
        num_projected = projected_mask.sum()
        if num_projected > 0:
            propagated_labels[seed_indices[projected_mask]] = seed_labels[projected_mask]
            print(f"Projected {num_projected} seed labels onto full cloud (tolerance: {seed_tolerance}m)")
        else:
            # Last resort: use a larger tolerance (up to 0.1m)
            print(f"Warning: No seed points matched even with {seed_tolerance}m tolerance. "
                  f"Trying with 0.1m tolerance...")
            seed_tolerance = 0.1
            projected_mask = seed_distances <= seed_tolerance
            num_projected = projected_mask.sum()
            if num_projected > 0:
                propagated_labels[seed_indices[projected_mask]] = seed_labels[projected_mask]
                print(f"Projected {num_projected} seed labels onto full cloud (tolerance: {seed_tolerance}m)")
            else:
                print(f"Error: No seed points matched even with 0.1m tolerance. "
                      f"Min distance: {seed_distances.min():.4f}m, Max distance: {seed_distances.max():.4f}m, "
                      f"Mean distance: {seed_distances.mean():.4f}m")
                print("This suggests the seed points and full cloud may be from different coordinate systems.")
    
    initial_labeled = (propagated_labels >= 1).sum()
    print(f"Initialized with {initial_labeled} labeled points from seeds")
    
    # Iterative wave propagation
    iteration = 0
    while iteration < max_iterations:
        wave_start_time = time.time()
        
        # Find currently labeled points
        labeled_mask = propagated_labels >= 1
        if not np.any(labeled_mask):
            print(f"Wave {iteration + 1}: No labeled points remaining, stopping.")
            break
        
        labeled_points = points[labeled_mask]
        labeled_labels = propagated_labels[labeled_mask]
        
        # Find unlabeled points
        unlabeled_mask = propagated_labels < 1
        unlabeled_points = points[unlabeled_mask]
        unlabeled_indices = np.where(unlabeled_mask)[0]
        
        if len(unlabeled_points) == 0:
            print(f"Wave {iteration + 1}: All points are labeled, stopping.")
            break
        
        # Build KDTree on currently labeled points
        kdtree_start = time.time()
        kdtree = cKDTree(labeled_points)
        kdtree_time = time.time() - kdtree_start
        
        # Query unlabeled points: find labeled points within wave_distance
        # Query up to 10 nearest neighbors (or all if fewer than 10)
        k_neighbors = min(10, len(labeled_points))
        
        query_start = time.time()
        # Handle case where there's only one unlabeled point (query returns 1D array)
        if len(unlabeled_points) == 1:
            distances, indices = kdtree.query(
                unlabeled_points,
                k=k_neighbors,
                distance_upper_bound=wave_distance
            )
            # Reshape to 2D for consistency
            if distances.ndim == 0:
                distances = np.array([[distances]])
                indices = np.array([[indices]])
            elif distances.ndim == 1:
                distances = distances.reshape(1, -1)
                indices = indices.reshape(1, -1)
        else:
            distances, indices = kdtree.query(
                unlabeled_points,
                k=k_neighbors,
                distance_upper_bound=wave_distance
            )
        query_time = time.time() - query_start
        
        # Apply strategy to handle multiple labeled points
        strategy_start = time.time()
        if multiple_label_strategy == 'hybrid':
            new_labels = _apply_hybrid_strategy(
                distances, indices, labeled_labels, wave_distance, weight_tolerance
            )
        elif multiple_label_strategy == 'closest':
            # Use closest labeled point (vectorized)
            new_labels = np.full(len(unlabeled_points), -1, dtype=int)
            # Vectorized: check which points have closest neighbor within threshold
            within_threshold = distances[:, 0] <= wave_distance
            # Vectorized: assign labels for points within threshold
            new_labels[within_threshold] = labeled_labels[indices[within_threshold, 0]]
        elif multiple_label_strategy == 'majority':
            # Use majority label (parallel with shared memory + vectorization)
            # Strategy: Split points into batches, process batches in parallel using shared memory
            # Each batch uses fully vectorized NumPy operations
            new_labels = np.full(len(unlabeled_points), -1, dtype=int)
            
            # Check if we should use parallel processing
            # Use threading if we're in a nested multiprocessing context, multiprocessing otherwise
            in_mp_context = _is_in_multiprocessing_context()
            use_parallel = (n_jobs is not None and n_jobs > 1 and len(unlabeled_points) > 10000)
            
            if use_parallel:
                # Split points into batches
                batch_size = max(50000, len(unlabeled_points) // (n_jobs * 2))  # Ensure at least 2 batches per worker
                batches = []
                for i in range(0, len(unlabeled_points), batch_size):
                    end_idx = min(i + batch_size, len(unlabeled_points))
                    batches.append((i, end_idx))
                
                parallel_success = False
                
                if in_mp_context:
                    # We're in a nested multiprocessing context - use threading
                    # Threads share memory, so we can pass arrays directly (no shared_memory needed)
                    try:
                        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
                            # Create batch info tuples for threading (pass arrays directly)
                            batch_args = [(start_idx, end_idx, distances, indices, labeled_labels, wave_distance)
                                        for start_idx, end_idx in batches]
                            results = list(executor.map(_process_batch_majority_thread, batch_args))
                        
                        # Combine results
                        for start_idx, end_idx, batch_labels in results:
                            new_labels[start_idx:end_idx] = batch_labels
                        parallel_success = True
                    except Exception as e:
                        print(f"Warning: Threading failed, falling back to sequential: {e}")
                        parallel_success = False
                        use_parallel = False
                else:
                    # We're in the main process - use multiprocessing with shared memory
                    # Create shared memory for read-only data
                    shm_dist = shared_memory.SharedMemory(create=True, size=distances.nbytes)
                    shm_dist_arr = np.ndarray(distances.shape, dtype=distances.dtype, buffer=shm_dist.buf)
                    shm_dist_arr[:] = distances[:]  # Copy once
                    
                    shm_idx = shared_memory.SharedMemory(create=True, size=indices.nbytes)
                    shm_idx_arr = np.ndarray(indices.shape, dtype=indices.dtype, buffer=shm_idx.buf)
                    shm_idx_arr[:] = indices[:]  # Copy once
                    
                    shm_lbl = shared_memory.SharedMemory(create=True, size=labeled_labels.nbytes)
                    shm_lbl_arr = np.ndarray(labeled_labels.shape, dtype=labeled_labels.dtype, buffer=shm_lbl.buf)
                    shm_lbl_arr[:] = labeled_labels[:]  # Copy once
                    
                    try:
                        # Create batch info tuples for multiprocessing (use shared memory names)
                        batch_args = []
                        for start_idx, end_idx in batches:
                            batch_args.append((
                                start_idx, end_idx,
                                shm_dist.name, shm_idx.name, shm_lbl.name,
                                distances.shape, indices.shape, labeled_labels.shape,
                                wave_distance
                            ))
                        
                        # Process batches in parallel
                        with mp.Pool(n_jobs) as pool:
                            results = pool.map(_process_batch_majority_shared, batch_args)
                        
                        # Combine results
                        for start_idx, end_idx, batch_labels in results:
                            new_labels[start_idx:end_idx] = batch_labels
                        parallel_success = True
                    except (RuntimeError, AssertionError) as e:
                        # If we can't create a pool, fall back to sequential
                        if "daemonic" in str(e) or "not allowed to have children" in str(e):
                            print(f"Warning: Multiprocessing failed, falling back to sequential: {e}")
                            parallel_success = False
                            use_parallel = False
                        else:
                            raise
                    finally:
                        # Cleanup shared memory (always, whether parallel succeeded or failed)
                        shm_dist.close()
                        shm_dist.unlink()
                        shm_idx.close()
                        shm_idx.unlink()
                        shm_lbl.close()
                        shm_lbl.unlink()
            
            # Sequential processing (fallback for small datasets, n_jobs=None/1, or nested multiprocessing)
            if not use_parallel:
                # Sequential processing (fallback for small datasets or n_jobs=None/1)
                # Use the same vectorized approach but without parallelization
                valid_mask = distances <= wave_distance  # Shape: (N, k)
                valid_neighbor_counts = np.sum(valid_mask, axis=1)  # Shape: (N,)
                
                # Process points with 1 neighbor (vectorized - fast path)
                single_mask = valid_neighbor_counts == 1
                if np.any(single_mask):
                    first_valid_idx = np.argmax(valid_mask[single_mask], axis=1)
                    point_indices_single = np.where(single_mask)[0]
                    neighbor_indices_single = indices[single_mask, first_valid_idx]
                    
                    max_valid_idx = len(labeled_labels) - 1
                    valid_single_mask = neighbor_indices_single <= max_valid_idx
                    if np.any(valid_single_mask):
                        new_labels[point_indices_single[valid_single_mask]] = labeled_labels[neighbor_indices_single[valid_single_mask]]
                
                # Process points with 2+ neighbors (vectorized grouping)
                multiple_mask = valid_neighbor_counts >= 2
                if np.any(multiple_mask):
                    N = len(unlabeled_points)
                    k = distances.shape[1]
                    point_indices = np.repeat(np.arange(N), k)
                    valid_flat = valid_mask.ravel()
                    
                    if np.any(valid_flat):
                        neighbor_indices_flat = indices.ravel()[valid_flat]
                        valid_point_flat = point_indices[valid_flat]
                        
                        multiple_point_indices = np.where(multiple_mask)[0]
                        multiple_neighbors_point_mask = np.isin(valid_point_flat, multiple_point_indices)
                        neighbor_indices_flat = neighbor_indices_flat[multiple_neighbors_point_mask]
                        valid_point_flat = valid_point_flat[multiple_neighbors_point_mask]
                        
                        max_valid_idx = len(labeled_labels) - 1
                        valid_idx_mask = neighbor_indices_flat <= max_valid_idx
                        neighbor_indices_flat = neighbor_indices_flat[valid_idx_mask]
                        valid_point_flat = valid_point_flat[valid_idx_mask]
                        
                        if len(neighbor_indices_flat) > 0:
                            neighbor_labels_flat = labeled_labels[neighbor_indices_flat]
                            sort_idx = np.argsort(valid_point_flat)
                            sorted_points = valid_point_flat[sort_idx]
                            sorted_labels = neighbor_labels_flat[sort_idx]
                            
                            unique_points, group_starts = np.unique(sorted_points, return_index=True)
                            group_ends = np.append(group_starts[1:], len(sorted_points))
                            
                            for i, point_idx in enumerate(unique_points):
                                start = group_starts[i]
                                end = group_ends[i]
                                point_labels = sorted_labels[start:end]
                                
                                if len(point_labels) > 0:
                                    max_label = point_labels.max()
                                    if max_label >= 0:
                                        counts = np.bincount(point_labels, minlength=max_label + 1)
                                        majority_label = np.argmax(counts)
                                        new_labels[point_idx] = majority_label
        elif multiple_label_strategy == 'same_only':
            # Only label if all nearby points have same label (optimized)
            new_labels = np.full(len(unlabeled_points), -1, dtype=int)
            for i in range(len(unlabeled_points)):
                within_threshold = distances[i] <= wave_distance
                if np.any(within_threshold):
                    nearby_labels = labeled_labels[indices[i][within_threshold]]
                    # Faster check: if min == max, all labels are the same
                    if len(nearby_labels) > 0 and nearby_labels.min() == nearby_labels.max():
                        new_labels[i] = nearby_labels[0]
        elif multiple_label_strategy == 'weighted':
            # Use weighted by inverse distance (partially vectorized)
            new_labels = np.full(len(unlabeled_points), -1, dtype=int)
            for i in range(len(unlabeled_points)):
                within_threshold = distances[i] <= wave_distance
                if np.any(within_threshold):
                    nearby_labels = labeled_labels[indices[i][within_threshold]]
                    nearby_distances = distances[i][within_threshold]
                    # Weight by inverse distance (vectorized)
                    weights = 1.0 / (nearby_distances + 1e-6)
                    # Use np.add.at for vectorized weight accumulation
                    max_label = nearby_labels.max()
                    if max_label >= 0:
                        label_weights = np.zeros(max_label + 1, dtype=float)
                        np.add.at(label_weights, nearby_labels, weights)
                        best_label = np.argmax(label_weights)
                        new_labels[i] = best_label
        else:
            raise ValueError(f"Unknown multiple_label_strategy: {multiple_label_strategy}")
        strategy_time = time.time() - strategy_start
        
        # Update labels
        newly_labeled = new_labels >= 1
        num_newly_labeled = newly_labeled.sum()
        
        if num_newly_labeled == 0:
            print(f"Wave {iteration + 1}: No new points labeled, stopping.")
            break
        
        propagated_labels[unlabeled_indices[newly_labeled]] = new_labels[newly_labeled]
        
        iteration += 1
        total_labeled = (propagated_labels >= 1).sum()
        wave_time = time.time() - wave_start_time
        print(f"Wave {iteration}: Labeled {num_newly_labeled} new points. "
              f"Total labeled: {total_labeled} ({100*total_labeled/len(points):.2f}%) | "
              f"Time: {wave_time:.2f}s (KDTree build: {kdtree_time:.2f}s, "
              f"Query: {query_time:.2f}s, Strategy: {strategy_time:.2f}s)")
        
        # Stop if too few new points (convergence)
        if num_newly_labeled < min_new_points:
            print(f"Wave {iteration}: Stopping - only {num_newly_labeled} new points "
                  f"labeled (threshold: {min_new_points})")
            break
    
    # Final statistics
    num_labeled = (propagated_labels >= 1).sum()
    print(f"Iterative Wave Propagation complete after {iteration} waves.")
    print(f"Labeled {num_labeled} points out of {points.shape[0]} total "
          f"({100*num_labeled/points.shape[0]:.2f}%)")
    
    # If an output file is specified, save the labels
    if out_file is not None:
        write_start_time = time.time()
        if out_file.endswith('.ply'):
            # Prepare vertex data (optimized: vectorized, no Python loop)
            prep_start = time.time()
            vertex_dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('label', 'i4')]
            vertex_array = np.empty(len(points), dtype=vertex_dtype)
            vertex_array['x'] = points[:, 0].astype('f4')
            vertex_array['y'] = points[:, 1].astype('f4')
            vertex_array['z'] = points[:, 2].astype('f4')
            vertex_array['label'] = propagated_labels.astype('i4')
            prep_time = time.time() - prep_start
            
            # Write to file
            write_file_start = time.time()
            ply_element = PlyElement.describe(vertex_array, 'vertex')
            PlyData([ply_element], text=False).write(out_file)
            write_file_time = time.time() - write_file_start
            
            total_write_time = time.time() - write_start_time
            print(f"Labels saved to {out_file}")
            print(f"  Write time: prep={prep_time:.2f}s, write_file={write_file_time:.2f}s, total={total_write_time:.2f}s")
        elif out_file.endswith('.pts'):
            # Prepare data
            prep_start = time.time()
            data = np.column_stack((points, propagated_labels))
            prep_time = time.time() - prep_start
            
            # Write to file
            write_file_start = time.time()
            np.savetxt(out_file, data, fmt=['%.3f', '%.3f', '%.3f', '%d'], delimiter=' ')
            write_file_time = time.time() - write_file_start
            
            total_write_time = time.time() - write_start_time
            print(f"Labels saved to {out_file}")
            print(f"  Write time: prep={prep_time:.2f}s, write_file={write_file_time:.2f}s, total={total_write_time:.2f}s")
        else:
            print(f"Unsupported output file format: {out_file}. Labels not saved.")
    else:
        print("No output file specified, labels not saved.")
    
    return propagated_labels


def run_label_propagation(
    infile: str,
    labeled_file: str,
    method: str = 'region_growing',
    **kwargs
) -> np.ndarray:
    """Runs label propagation on a point cloud using the specified method.

    This high-level function loads the necessary point cloud files and dispatches
    to the appropriate labeling algorithm.

    Args:
        infile (str): Path to the input point cloud file to be labeled (.pts or .ply).
        labeled_file (str): Path to the file containing the labeled seed points
            (.pts or .ply with X,Y,Z,label).
        method (str): The method to use for label propagation. One of
            'region_growing', 'region_growing_layered', 'layered_nn', 'distance_based', or 'iterative_distance_based'. 
            Defaults to 'region_growing'.
            - 'region_growing': Sequential Z-ordered region growing (original)
            - 'region_growing_layered': Z-layer batched region growing (faster, maintains Z-ordering)
            - 'layered_nn': Layered nearest neighbor method
            - 'distance_based': Distance-based nearest neighbor (fastest for dense clouds)
            - 'iterative_distance_based': Iterative wave propagation with distance-based queries
        **kwargs: Additional keyword arguments to be passed to the chosen
            labeling function (e.g., `search_radius`, `layer_height`).

    Returns:
        np.ndarray: The array of labels for the input point cloud.

    Raises:
        ValueError: If an unsupported file format or method is provided.
        FileNotFoundError: If either input file cannot be found.
    """
    # Load the point clouds
    if infile.endswith('.ply'):
        ply_data = PlyData.read(infile)
        points = np.vstack([ply_data['vertex'][col] for col in ['x', 'y', 'z']]).T
    elif infile.endswith('.pts'):
        points = np.loadtxt(infile)
        # If points have more than 3 columns, keep only XYZ
        if points.shape[1] > 3:
            points = points[:, :3]
    elif infile.endswith('.las') or infile.endswith('.laz'):
        las = laspy.read(infile)
        points = np.vstack((las.x, las.y, las.z)).T
    else:
        raise ValueError(f"Unsupported input file format: {infile}")
    
    print(f"Loaded {points.shape[0]} points from {infile}")
    
    # Load the labeled seed points
    if labeled_file.endswith('.ply'):
        ply_data = PlyData.read(labeled_file)
        labeled_seed_points = np.vstack([ply_data['vertex'][col] for col in ['x', 'y', 'z', 'label']]).T
    elif labeled_file.endswith('.pts'):
        labeled_seed_points = np.loadtxt(labeled_file)
    else:
        raise ValueError(f"Unsupported labeled seed points file format: {labeled_file}")
    
    print(f"Loaded {labeled_seed_points.shape[0]} labeled seed points from {labeled_file}")
    
    if method == 'region_growing':
        return label_points_region_growing(points, labeled_seed_points, **kwargs)
    elif method == 'region_growing_layered':
        return label_points_region_growing_layered(points, labeled_seed_points, **kwargs)
    elif method == 'layered_nn':
        return label_points_layered_nn(points, labeled_seed_points, **kwargs)
    elif method == 'distance_based':
        return label_points_distance_based(points, labeled_seed_points, **kwargs)
    elif method == 'iterative_distance_based':
        return label_points_iterative_distance_based(points, labeled_seed_points, **kwargs)
    else:
        raise ValueError(f"Unknown method: {method}. Choose from 'region_growing', 'region_growing_layered', 'layered_nn', 'distance_based', or 'iterative_distance_based'")

# ----- old code -----
# Not suggested.
from collections import defaultdict
def label_pts_from_core_example(pts, edges, lpts):
    # pts = defaultdict(lambda: -1)  # (x,y,z):id
    # edges = defaultdict(lambda: set())  # pid1: set of ALL neighbors pid2. pid2 can < pid
    # edges can be made from the pts via alpha shape
    # lpts: list of (x,y,z,l)

    # todo: we can add the confidence for each point here.
    # lpts confidence = 1
    # new labeled points = 0.5, or decreased based on the distance to nearest 1-confidence & same label pts
    pts_num = len(pts)
    out_pts = [None for i in range(pts_num)]
    for p in pts:
        x, y, z = p
        pid = pts[p]
        out_pts[pid] = (x, y, z)
    sort_pts = sorted(out_pts, key=lambda x: x[2])
    lbl_sorted_pids = [-1 for i in range(pts_num)]
    old2new = defaultdict(lambda: -1)
    for i in range(pts_num):
        p = sort_pts[i]
        oldpid = pts[p]
        old2new[oldpid] = i
    # init
    for lp in lpts:
        x, y, z, l = lp
        p = (x, y, z)
        if p in pts:
            oldpid = pts[p]
            newpid = old2new[oldpid]
            lbl_sorted_pids[newpid] = l
    # iterate
    labelnew = True
    while labelnew:
        labelnew = False
        for i in range(pts_num):
            l = lbl_sorted_pids[i]
            if l > -1:
                p = sort_pts[i]
                pid1 = pts[p]
                for np in edges[pid1]:
                    nnp = old2new[np]
                    if lbl_sorted_pids[nnp] < 0:
                        lbl_sorted_pids[nnp] = l
                        labelnew = True

    #
    final_pts = list()
    for i in range(pts_num):
        x, y, z = sort_pts[i]
        l = lbl_sorted_pids[i]
        final_pts.append((x, y, z, l))
    return final_pts