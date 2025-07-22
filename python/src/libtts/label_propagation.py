import numpy as np
from scipy.spatial import KDTree

import heapq

def label_points_layered_nn(points, labeled_seed_points, layer_height=1.0, max_search_radius=1.0, out_file=None):
    """
    Labels unlabeled points from labeled points using a layered Nearest Neighbor approach with a distance cap.
    Args:
        points (np.ndarray): The target Nx3 point cloud to be labeled.
        labeled_seed_points (np.ndarray): An Mx4 array of initially labeled points (x, y, z, label).
        layer_height (float): The thickness of each horizontal layer.
        max_search_radius (float): The maximum distance to a labeled point for propagation.

    Returns:
        np.ndarray: The array of labels for the target 'points' cloud.
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
        np.savetxt(out_file, np.column_stack((points, propagated_labels)), fmt='%.3f', delimiter=' ')
        print(f"Labels saved to {out_file}")
    return propagated_labels


def label_points_region_growing(points, labeled_seed_points, search_radius=0.5, out_file=None):
    """
    Labels unlabeled points using a region growing method ordered by Z-height.
    This version is robust to mismatches between the main point cloud and seed points.

    Args:
        points (np.ndarray): The target Nx3 point cloud to be labeled.
        labeled_seed_points (np.ndarray): An Mx4 array of initially labeled points (x, y, z, label).
        search_radius (float): The distance to consider points as "connected".
        out_file (str, optional): If provided, saves the final labeled points to this path.

    Returns:
        np.ndarray: The array of labels for the target 'points' cloud.
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
    
    # todo: support .ply 
    
    # If an output file is specified, save the labels: xyzl
    if out_file is not None:
        np.savetxt(out_file, np.column_stack((points, propagated_labels)), fmt=['%.3f', '%.3f', '%.3f', '%d'], delimiter=' ')
        print(f"Labels saved to {out_file}")
        
    return propagated_labels

# ----- 
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