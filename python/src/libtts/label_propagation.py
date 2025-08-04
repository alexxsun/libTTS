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

import numpy as np
from scipy.spatial import KDTree
import heapq
from plyfile import PlyData, PlyElement

def label_points_layered_nn(points, labeled_seed_points, layer_height=1.0, max_search_radius=1.0, out_file=None):
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
        np.ndarray: An array of shape (N,) containing the propagated labels for
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
            vertices = np.column_stack((points, propagated_labels))
            vertex_element = PlyElement.describe(vertices, 'vertex', comments=['Labeled points'])
            ply_data = PlyData([vertex_element])
            ply_data.write(out_file)
            print(f"Labels saved to {out_file}")
        elif out_file.endswith('.pts'):
            np.savetxt(out_file, np.column_stack((points, propagated_labels)), fmt='%.3f', delimiter=' ')
            print(f"Labels saved to {out_file}")
        else:
            print(f"Unsupported output file format: {out_file}. Labels not saved.")
    return propagated_labels


def label_points_region_growing(points, labeled_seed_points, search_radius=0.5, out_file=None):
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
        np.ndarray: An array of shape (N,) containing the propagated labels for
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
            vertices = np.column_stack((points, propagated_labels))
            vertex_element = PlyElement.describe(vertices, 'vertex', comments=['Labeled points'])
            ply_data = PlyData([vertex_element])
            ply_data.write(out_file)
            print(f"Labels saved to {out_file}")
        elif out_file.endswith('.pts'):
            np.savetxt(out_file, np.column_stack((points, propagated_labels)), fmt=['%.3f', '%.3f', '%.3f', '%d'], delimiter=' ')
            print(f"Labels saved to {out_file}")
        else:
            print(f"Unsupported output file format: {out_file}. Labels not saved.")
    else:
        print("No output file specified, labels not saved.")
        
    return propagated_labels

def run_label_propagation(infile, labeled_file, method='region_growing', **kwargs):
    """Runs label propagation on a point cloud using the specified method.

    This high-level function loads the necessary point cloud files and dispatches
    to the appropriate labeling algorithm.

    Args:
        infile (str): Path to the input point cloud file to be labeled (.pts or .ply).
        labeled_file (str): Path to the file containing the labeled seed points
            (.pts or .ply with X,Y,Z,label).
        method (str): The method to use for label propagation. One of
            'region_growing' or 'layered_nn'. Defaults to 'region_growing'.
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
    else:
        raise ValueError(f"Unsupported input file format: {infile}")
    
    # Load the labeled seed points
    if labeled_file.endswith('.ply'):
        ply_data = PlyData.read(labeled_file)
        labeled_seed_points = np.vstack([ply_data['vertex'][col] for col in ['x', 'y', 'z', 'label']]).T
    elif labeled_file.endswith('.pts'):
        labeled_seed_points = np.loadtxt(labeled_file)
    else:
        raise ValueError(f"Unsupported labeled seed points file format: {labeled_file}")
    
    if method == 'region_growing':
        return label_points_region_growing(points, labeled_seed_points, **kwargs)
    elif method == 'layered_nn':
        return label_points_layered_nn(points, labeled_seed_points, **kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")

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