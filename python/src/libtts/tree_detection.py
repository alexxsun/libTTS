"""Functions for detecting potential tree locations from a normalized point cloud.

Example:
# simple tree detection
tree_locfile = libtts.run_tree_detection(infile = tls_veg_file, outfile= "tree_locations.pts", 
                                         method='base',
                                         height_min=0.5, height_max=1.0, 
                                         max_avg_dist=1.0,
                                         eps=0.2, min_samples=5)

# geometry-based tree detection
tree_locfile = libtts.run_tree_detection(infile = out_ply_file, outfile= "tree_locations.pts", 
                                         method='geometry',
                                         height_min=0.5, height_max=1.0, 
                                         n_neighbors_pca = 20,
                                         max_linearity=0.2, max_knn_dist = 0.02,
                                         eps=0.1, min_samples=20)

# gridding-based tree detection
tree_locfile = libtts.run_tree_detection(infile = ptsfile, outfile= "tree_locations_gridding.pts", 
                                         method='grid',
                                         height_min=0.5, height_max=1.0, 
                                         grid_size = 0.05,
                                         eps=0.1, min_samples=2,
                                         does_plot=True)

"""
import argparse
import os
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from typing import Optional, Union

try:
    from plyfile import PlyData, PlyElement
    PLYFILE_ENABLED = True
except ImportError:
    PLYFILE_ENABLED = False

try:
    import matplotlib.pyplot as plt
    PLOTTING_ENABLED = True
except ImportError:
    PLOTTING_ENABLED = False


# --- Utility Functions ---

def _load_points(pts_file: str) -> np.ndarray:
    """Loads points from a .pts, .txt, or .ply file into an NxM numpy array."""
    if not os.path.exists(pts_file):
        raise FileNotFoundError(f"Input file not found: {pts_file}")

    if pts_file.lower().endswith((".pts", ".txt")):
        points = np.loadtxt(pts_file)
    elif pts_file.lower().endswith(".ply"):
        if not PLYFILE_ENABLED:
            raise ImportError("The 'plyfile' library is required to read .ply files.")
        ply_data = PlyData.read(pts_file)
        vertices = ply_data['vertex']
        data = [vertices['x'], vertices['y'], vertices['z'], vertices['h']]
        points = np.vstack(data).T
    else:
        raise ValueError("Unsupported file format. Please use .pts, .txt, or .ply.")

    if points.shape[1] == 3:
        print(f"xyz -> xyzz: repeating z as h.")
        points = np.column_stack((points, points[:, 2]))
    return points

def _save_points(points: np.ndarray, outfile: str):
    """Saves points to a file, detecting format from extension."""
    if outfile.lower().endswith((".pts", ".txt")):
        np.savetxt(outfile, points, fmt="%.3f")
    elif outfile.lower().endswith(".ply"):
        if not PLYFILE_ENABLED:
            raise ImportError("The 'plyfile' library is required to save .ply files.")
        dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]
        if points.shape[1] > 3:
            dtype.append(('label', 'i4'))
        structured_points = np.array([tuple(row) for row in points], dtype=dtype)
        el = PlyElement.describe(structured_points, 'vertex')
        PlyData([el]).write(outfile)
    else:
        raise ValueError("Unsupported output format. Please use .pts, .txt, or .ply.")
    print(f"Saved {len(points)} points to {outfile}")

# --- Core Algorithm Functions ---

def filter_tree_bases(
    points_xyzh: np.ndarray,
    height_min: float = 0.5,
    height_max: float = 1.0,
    knn: int = 4,
    max_avg_dist: float = 0.1
) -> np.ndarray:
    """
    Filters a point cloud to find dense points in a low height range.
    """
    heights = points_xyzh[:, 3]
    slice_mask = (heights >= height_min) & (heights <= height_max)
    low_points = points_xyzh[slice_mask]

    if low_points.shape[0] <= knn:
        print("Warning: Not enough points in the height slice for KNN filtering.")
        return np.array([])

    search_coords = low_points[:, [0, 1, 3]]
    nbrs = NearestNeighbors(n_neighbors=knn + 1, algorithm='kd_tree').fit(search_coords)
    distances, _ = nbrs.kneighbors(search_coords)
    avg_knn_dist = np.mean(distances[:, 1:], axis=1)
    dense_mask = avg_knn_dist < max_avg_dist
    filtered_points = low_points[dense_mask]
    
    print(f"Filtered {len(low_points)} slice points down to {len(filtered_points)} dense points.")
    return filtered_points


def filter_points_by_geometry(
    points_xyzh: np.ndarray,
    height_min: float = 0.5,
    height_max: float = 1.5,
    n_neighbors_pca: int = 15,
    max_linearity: float = 0.7,
    max_knn_dist: float = 0.1,
    knn_for_dist: int = 4
) -> np.ndarray:
    """
    Filters points based on local geometric properties (density and linearity)
    using a memory-efficient k-Nearest Neighbors approach.

    Args:
        points_xyzh: Input Nx4 point cloud.
        height_min: Minimum height for the analysis slice.
        height_max: Maximum height for the analysis slice.
        n_neighbors_pca: The number of neighbors to define the local neighborhood for PCA.
        max_linearity: Threshold for linearity (1D structure). Points below this are kept.
        max_knn_dist: Threshold for local density. Points below this are kept.
        knn_for_dist: Number of neighbors for the final density check.

    Returns:
        Filtered points that are dense and not strongly linear.
    """
    # 1. Select points within the specified height slice
    heights = points_xyzh[:, 3]
    slice_mask = (heights >= height_min) & (heights <= height_max)
    low_points = points_xyzh[slice_mask]

    if low_points.shape[0] < n_neighbors_pca:
        print("Warning: Not enough points in height slice for geometry filtering.")
        return np.array([])
    

    # debug: save low_points for inspection
    debug_dir = "/home/alex/Projects/libTTS_public/some_examples/"
    np.savetxt(f"{debug_dir}/debug_low_points.pts", low_points, fmt="%.3f")
    print(f"# pts: {len(low_points)} in height slice {height_min} - {height_max}.")

    # 2. Find neighbors for all points at once
    search_coords = low_points[:, :3] # Use XYZ for geometry
    nbrs = NearestNeighbors(n_neighbors=n_neighbors_pca).fit(search_coords)
    distances, indices = nbrs.kneighbors(search_coords)

    # 3. Iterate and calculate geometric features for each point's neighborhood
    points_to_keep_mask = np.zeros(low_points.shape[0], dtype=bool)
    # debug: save all linearity and density values for inspection
    linearity_values = []
    density_values = []
    for i in range(low_points.shape[0]):
        neighbor_indices = indices[i]        
        neighborhood_points = search_coords[neighbor_indices]
        
        # a) Calculate linearity via PCA
        pca = PCA(n_components=3)
        pca.fit(neighborhood_points)
        eigenvalues = pca.explained_variance_
        
        # Linearity is (e1 - e2) / e1, where e1 is the largest eigenvalue
        linearity = (eigenvalues[0] - eigenvalues[1]) / eigenvalues[0] if eigenvalues[0] > 0 else 0

        # b) Calculate local density (KNN distance)
        avg_knn_dist = np.mean(distances[i, 1:knn_for_dist+1])

        #
        linearity_values.append(linearity)
        density_values.append(avg_knn_dist)

        # c) Apply filters
        if linearity < max_linearity and avg_knn_dist < max_knn_dist:
            points_to_keep_mask[i] = True
            
    # debug: save x,y,z, linearity and density values for inspection
    linearity_values = np.array(linearity_values)
    density_values = np.array(density_values)
    debug_data = np.column_stack((low_points[:, :3], linearity_values, density_values))
    np.savetxt(f"{debug_dir}/debug_geometry_values.pts", debug_data, fmt="%.3f")
    
    filtered_points = low_points[points_to_keep_mask]
    print(f"Filtered {len(low_points)} slice points down to {len(filtered_points)} based on geometry.")
    return filtered_points


def cluster_points_dbscan(
    points: np.ndarray,
    eps: float = 0.2,
    min_samples: int = 10,
    use_2d: bool = True
) -> np.ndarray:
    """
    Clusters points using DBSCAN and returns points with cluster labels.
    """
    if points.shape[0] == 0:
        return np.array([])

    coords_to_cluster = points[:, :2] if use_2d else points[:, :3]
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords_to_cluster)
    labels = db.labels_
    non_noise_mask = labels != -1
    
    if not np.any(non_noise_mask):
        return np.array([])

    labeled_points = points[non_noise_mask]
    final_labels = labels[non_noise_mask] + 1
    
    print(f"Found {len(np.unique(final_labels))} clusters.")
    return np.column_stack((labeled_points, final_labels))


# --- High-Level API Functions ---

def detect_trees_by_base_density(points_xyzh: np.ndarray, **kwargs) -> np.ndarray:
    """Workflow: Filters for dense points at a low height and clusters them."""
    filter_params = {
        "height_min": kwargs.get("height_min", 0.5),
        "height_max": kwargs.get("height_max", 1.5),
        "knn": kwargs.get("knn", 4),
        "max_avg_dist": kwargs.get("max_avg_dist", 0.1)
    }
    cluster_params = {
        "eps": kwargs.get("eps", 0.2),
        "min_samples": kwargs.get("min_samples", 10)
    }
    tree_bases = filter_tree_bases(points_xyzh, **filter_params)
    return cluster_points_dbscan(tree_bases, **cluster_params)


def detect_trees_by_geometry(points_xyzh: np.ndarray, **kwargs) -> np.ndarray:
    """Workflow: Filters points by local geometry (linearity/density) and clusters them."""
    filter_params = {
        "height_min": kwargs.get("height_min", 0.5),
        "height_max": kwargs.get("height_max", 1.5),
        "n_neighbors_pca": kwargs.get("n_neighbors_pca", 15),
        "max_linearity": kwargs.get("max_linearity", 0.7),
        "max_knn_dist": kwargs.get("max_knn_dist", 0.1)
    }
    cluster_params = {
        "eps": kwargs.get("eps", 0.2),
        "min_samples": kwargs.get("min_samples", 10)
    }
    seed_points = filter_points_by_geometry(points_xyzh, **filter_params)
    return cluster_points_dbscan(seed_points, **cluster_params)

def detect_trees_by_gridding(points_xyzh: np.ndarray, **kwargs) -> np.ndarray:
    """Workflow: Creates a 2D histogram, filters dense cells, and clusters them."""
    # 1. Extract parameters
    height_min = kwargs.get("height_min", 0.5)
    height_max = kwargs.get("height_max", 1.0)
    grid_size = kwargs.get("grid_size", 0.05)
    min_points_per_cell = kwargs.get("min_points_per_cell", 100)
    eps = kwargs.get("eps", 0.1)
    min_samples = kwargs.get("min_samples", 2)
    does_plot = kwargs.get("does_plot", False)
    if does_plot and not PLOTTING_ENABLED:
        raise ImportError("The 'matplotlib' library is required for plotting.")

    # 2. Select points within the specified height slice
    heights = points_xyzh[:, 3]
    slice_mask = (heights >= height_min) & (heights <= height_max)
    low_points = points_xyzh[slice_mask]

    if low_points.shape[0] == 0:
        print("Warning: No points in the height slice for gridding.")
        return np.array([])

    # 3. Create 2D histogram
    x_min, x_max = low_points[:, 0].min(), low_points[:, 0].max()
    y_min, y_max = low_points[:, 1].min(), low_points[:, 1].max()
    x_bins = np.arange(x_min, x_max + grid_size, grid_size)
    y_bins = np.arange(y_min, y_max + grid_size, grid_size)
    hist, x_edges, y_edges = np.histogram2d(low_points[:, 0], low_points[:, 1], bins=[x_bins, y_bins])

    # 4. Filter for dense cells
    dense_cell_indices = np.argwhere(hist >= min_points_per_cell)
    if dense_cell_indices.shape[0] == 0:
        print("No dense cells found above the threshold.")
        return np.array([])
    
    if does_plot:
        # Plot the histogram with dense cells highlighted
        plt.figure(figsize=(8, 6))
        plt.imshow(hist.T, origin='lower', cmap='viridis', extent=[x_min, x_max, y_min, y_max], aspect='auto')
        for cell in dense_cell_indices:
            plt.gca().add_patch(plt.Rectangle((x_edges[cell[0]], y_edges[cell[1]]), grid_size, grid_size, edgecolor='red', facecolor='none', lw=1, alpha = 0.9))
        plt.colorbar(label='Number of Points')
        plt.title('Selected Grid Cells with More than {} Points'.format(min_points_per_cell))
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        #plt.savefig('selected_grid_cells.png')
        plt.show()
        plt.close()
        
        # show statistics of the histogram
        hist2 = hist[hist > 50]  # remove small cells
        print("Histogram statistics: >50 points per cell")
        print(f"Shape: {hist2.shape}")
        print(f"Total points: {hist2.sum()}")
        print(f"Max points in a cell: {hist2.max()}")
        print(f"Min points in a cell: {hist2.min()}")
        print(f"Mean points in a cell: {hist2.mean()}")
        print(f"Std points in a cell: {hist2.std()}")

    # 5. Cluster the dense cells
    # Get the XY coordinates of the bottom-left corner of each dense cell
    dense_cell_coords = np.array([
        [x_edges[idx[0]], y_edges[idx[1]]] for idx in dense_cell_indices
    ])
    
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(dense_cell_coords)
    labels = db.labels_
    
    # 6. Label original points based on the cell clusters
    # Create a mapping from cell index to cluster label
    cell_idx_to_label = {tuple(cell_idx): label + 1 for cell_idx, label in zip(dense_cell_indices, labels) if label != -1}

    if not cell_idx_to_label:
        print("No valid clusters found after DBSCAN.")
        return np.array([])

    # Vectorized approach to find which cell each point belongs to
    point_x_indices = ((low_points[:, 0] - x_min) / grid_size).astype(int)
    point_y_indices = ((low_points[:, 1] - y_min) / grid_size).astype(int)
    point_cell_indices = np.vstack([point_x_indices, point_y_indices]).T
    
    # Create a final label array
    final_labels = np.full(low_points.shape[0], -1, dtype=int)
    for i in range(low_points.shape[0]):
        cell_tuple = tuple(point_cell_indices[i])
        if cell_tuple in cell_idx_to_label:
            final_labels[i] = cell_idx_to_label[cell_tuple]
            
    # Filter for points that received a valid label
    valid_mask = final_labels != -1
    labeled_points = low_points[valid_mask]
    final_labels = final_labels[valid_mask]
    
    print(f"Found {len(np.unique(final_labels))} clusters using the gridding method.")
    return np.column_stack((labeled_points, final_labels))

def run_tree_detection(
    infile: str,
    outfile: Optional[str] = None,
    method: str = 'base',
    **kwargs
) -> Optional[Union[np.ndarray, str]]:
    """
    Runs a complete tree detection workflow on a file.

    If an outfile is provided, it saves the result and returns the output path.
    If outfile is None, it returns the labeled points as a NumPy array.

    Args:
        infile (str): Path to the input point cloud file.
        outfile (str, optional): Path to save the output labeled tree points. Defaults to None.
        method (str): The detection method to use, either 'base' or 'geometry'.
        **kwargs: Keyword arguments for the underlying detection functions.

    Returns:
        Optional[Union[np.ndarray, str]]: A NumPy array of labeled points (x,y,z,label)
                                           if outfile is None, or the output file path if
                                           outfile is provided. Returns None on failure.
    """
    try:
        points = _load_points(infile)
    except (FileNotFoundError, ValueError, ImportError) as e:
        print(f"Error loading file: {e}")
        return None

    if method == 'base':
        labeled_trees = detect_trees_by_base_density(points, **kwargs)
    elif method == 'geometry':
        labeled_trees = detect_trees_by_geometry(points, **kwargs)
    elif method == 'grid':
        labeled_trees = detect_trees_by_gridding(points, **kwargs)
    else:
        raise ValueError(f"Unknown detection method: '{method}'. Choose 'base', 'geometry', or 'grid'.")

    if labeled_trees.shape[0] == 0:
        print("No trees were detected.")
        return np.array([]) if outfile is None else None

    # Prepare data for output (x, y, z, label)
    output_data = labeled_trees[:, [0, 1, 2, 4]]

    if outfile:
        _save_points(output_data, outfile)
        return outfile
    else:
        return output_data


# --- Main CLI ---

def main():
    """Main function to run tree detection from the command line."""
    parser = argparse.ArgumentParser(description="Detects tree locations from a normalized point cloud.")
    parser.add_argument("-i", "--infile", required=True, help="Input point cloud file (.pts, .txt, .ply).")
    parser.add_argument("-o", "--outfile", required=True, help="Path to save the output labeled tree points.")
    parser.add_argument("--method", choices=['base', 'geometry', 'grid'], default='base', help="Detection method to use.")
    # Add other parameters...
    args = parser.parse_args()

    result = run_tree_detection(
        infile=args.infile,
        outfile=args.outfile,
        method=args.method
    )

    if result is not None and os.path.exists(str(result)):
        print(f"Tree detection complete. Output saved to: {result}")
    elif result is not None:
        print(f"Tree detection complete. Returned {len(result)} points.")
    else:
        print("Tree detection failed or produced no output.")

if __name__ == "__main__":
    main()