"""
Functions for detecting potential tree locations from a normalized point cloud.
Current implmentation: filtering and clustering points at a low height above the ground.
"""
import argparse
import os
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN

try:
    from plyfile import PlyData, PlyElement
    PLYFILE_ENABLED = True
except ImportError:
    PLYFILE_ENABLED = False

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
        # Dynamically build the array based on available fields (xyz are required)
        data = [vertices['x'], vertices['y'], vertices['z']]
        if 'h' in vertices.fields:
            data.append(vertices['h'])
        points = np.vstack(data).T
    else:
        raise ValueError("Unsupported file format. Please use .pts, .txt, or .ply.")

    # If only xyz is present, duplicate z as h for consistency
    if points.shape[1] == 3:
        points = np.column_stack((points, points[:, 2]))

    return points

def _save_points(points: np.ndarray, outfile: str):
    """Saves points to a file, detecting format from extension."""
    if outfile.lower().endswith((".pts", ".txt")):
        np.savetxt(outfile, points, fmt="%.3f")
    elif outfile.lower().endswith(".ply"):
        if not PLYFILE_ENABLED:
            raise ImportError("The 'plyfile' library is required to save .ply files.")
        # Create a structured array for plyfile
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
    Filters a point cloud to find dense points in a low height range,
    which likely correspond to tree bases.

    Args:
        points_xyzh (np.ndarray): Nx4 array of (x, y, z, h) points.
        height_min (float): The minimum height above ground to consider.
        height_max (float): The maximum height above ground to consider.
        knn (int): The number of nearest neighbors to use for denoising.
        max_avg_dist (float): The maximum average distance to neighbors for a
                              point to be kept.

    Returns:
        np.ndarray: A filtered array of points representing potential tree bases.
    """
    # 1. Select points within the specified height slice
    heights = points_xyzh[:, 3]
    slice_mask = (heights >= height_min) & (heights <= height_max)
    low_points = points_xyzh[slice_mask]

    if low_points.shape[0] <= knn:
        print("Warning: Not enough points in the height slice to perform filtering.")
        return np.array([])

    # 2. Denoise using KNN distance (vectorized for performance)
    # We use xyh for distance calculation as it's more robust than xyz
    search_coords = low_points[:, [0, 1, 3]]
    nbrs = NearestNeighbors(n_neighbors=knn + 1, algorithm='kd_tree').fit(search_coords)
    distances, _ = nbrs.kneighbors(search_coords)

    # Calculate the average distance to the k nearest neighbors (excluding the point itself)
    avg_knn_dist = np.mean(distances[:, 1:], axis=1)

    # Keep only the points that are in dense regions
    dense_mask = avg_knn_dist < max_avg_dist
    filtered_points = low_points[dense_mask]
    
    print(f"Filtered {len(low_points)} slice points down to {len(filtered_points)} dense points.")
    return filtered_points


def cluster_points_dbscan(
    points: np.ndarray,
    eps: float = 0.2,
    min_samples: int = 10,
    use_2d: bool = True
) -> np.ndarray:
    """
    Clusters points using DBSCAN and returns points with cluster labels.

    Args:
        points (np.ndarray): The NxM input points to cluster.
        eps (float): The maximum distance between two samples for one to be
                     considered as in the neighborhood of the other.
        min_samples (int): The number of samples in a neighborhood for a point
                           to be considered as a core point.
        use_2d (bool): If True, performs clustering on XY coordinates only.
                       Otherwise, uses XYZ.

    Returns:
        np.ndarray: An array of points with a new label column (x, y, z, ..., label).
                    Noise points are excluded. Labels start from 1.
    """
    if points.shape[0] == 0:
        return np.array([])

    coords_to_cluster = points[:, :2] if use_2d else points[:, :3]
    
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords_to_cluster)
    labels = db.labels_

    # Filter out noise points (label == -1)
    non_noise_mask = labels != -1
    labeled_points = points[non_noise_mask]
    final_labels = labels[non_noise_mask]

    # Increment labels so they start from 1 (a common convention)
    final_labels += 1
    
    print(f"Found {len(np.unique(final_labels))} clusters.")
    return np.column_stack((labeled_points, final_labels))


# --- High-Level API Function ---

def detect_trees(
    points_xyzh: np.ndarray,
    **kwargs
) -> np.ndarray:
    """
    Detects tree locations from a normalized point cloud.

    This workflow filters for dense points at a low height (tree bases)
    and then clusters them to identify individual trees.

    Args:
        points_xyzh (np.ndarray): The input Nx4 point cloud (x, y, z, h).
        **kwargs: Keyword arguments for the underlying functions:
            - height_min, height_max, knn, max_avg_dist (for filter_tree_bases)
            - eps, min_samples, use_2d (for cluster_points_dbscan)

    Returns:
        np.ndarray: An Mx5 array of labeled tree base points (x, y, z, h, label).
    """
    # Extract parameters for each step, providing defaults
    filter_params = {
        "height_min": kwargs.get("height_min", 0.5),
        "height_max": kwargs.get("height_max", 1.5),
        "knn": kwargs.get("knn", 4),
        "max_avg_dist": kwargs.get("max_avg_dist", 0.1)
    }
    cluster_params = {
        "eps": kwargs.get("eps", 0.2),
        "min_samples": kwargs.get("min_samples", 10),
        "use_2d": kwargs.get("use_2d", True)
    }

    tree_bases = filter_tree_bases(points_xyzh, **filter_params)
    labeled_trees = cluster_points_dbscan(tree_bases, **cluster_params)
    
    return labeled_trees

# --- Complete Workflow ---
from typing import Optional, Union

def run_tree_detection(
    infile: str,
    outfile: Optional[str] = None,
    height_min: float = 0.5,
    height_max: float = 1.5,
    max_dist: float = 0.1,
    eps: float = 0.2,
    min_samples: int = 10
)  -> Optional[np.ndarray]:
    """
    Runs the complete tree detection workflow from input file to output file.

    Args:
        infile (str): Path to the input point cloud file.
        outfile (str): Path to save the output labeled tree points.
        height_min (float): Minimum height of the slice to analyze.
        height_max (float): Maximum height of the slice to analyze.
        max_dist (float): Maximum average KNN distance for denoising.
        eps (float): DBSCAN eps parameter.
        min_samples (int): DBSCAN min_samples parameter.
    Returns:
        Optional[np.ndarray]: A NumPy array of labeled points (x,y,z,label) if outfile is None,
                              otherwise None.
    """
    points = _load_points(infile)

    if points.shape[1] == 3: # ensure points have height
        print("Input points only have x, y, z. Duplicating z as h for consistency.")
        points = np.column_stack((points, points[:, 2]))

    labeled_trees = detect_trees(
        points,
        height_min=height_min,
        height_max=height_max,
        max_avg_dist=max_dist,
        eps=eps,
        min_samples=min_samples
    )

    if labeled_trees.shape[0] > 0:
        # Conditionally save the file or return the data
        if outfile:
            # File-processing mode
            _save_points(labeled_trees[:, [0, 1, 2, 4]], outfile)
            return outfile # Return the output file path
        else:
            # Data-in, Data-out mode
            return labeled_trees[:, [0, 1, 2, 4]]
    else:
        print("No trees were detected.")

# --- Main CLI ---

def main():
    """Main function to run tree detection from the command line."""
    parser = argparse.ArgumentParser(description="Detects tree locations from a normalized point cloud.")
    parser.add_argument("infile", help="Path to the input point cloud file (.pts, .txt, or .ply with x,y,z,h fields).")
    parser.add_argument("outfile", help="Path to save the output labeled tree points.")
    parser.add_argument("--height_min", type=float, default=0.5, help="Minimum height of the slice to analyze.")
    parser.add_argument("--height_max", type=float, default=1.5, help="Maximum height of the slice to analyze.")
    parser.add_argument("--max_dist", type=float, default=0.1, help="Maximum average KNN distance for denoising.")
    parser.add_argument("--eps", type=float, default=0.2, help="DBSCAN eps parameter.")
    parser.add_argument("--min_samples", type=int, default=10, help="DBSCAN min_samples parameter.")

    args = parser.parse_args()

    try:
        points = _load_points(args.infile)
    except (FileNotFoundError, ValueError, ImportError) as e:
        print(f"Error: {e}")
        return

    # Run the main detection workflow
    labeled_trees = detect_trees(
        points,
        height_min=args.height_min,
        height_max=args.height_max,
        max_avg_dist=args.max_dist,
        eps=args.eps,
        min_samples=args.min_samples
    )

    if labeled_trees.shape[0] > 0:
        # Save the xyzl output
        _save_points(labeled_trees[:, [0, 1, 2, 4]], args.outfile)
    else:
        print("No trees were detected.")

    print("Tree detection complete.")

if __name__ == "__main__":
    main()
