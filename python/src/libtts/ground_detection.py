"""Functions for detecting ground points and classifying points (ground and vegetation).

Current implementation: creating a Digital Terrain Model (DTM), and classifying points as ground or vegetation based on the height.

Example:
    out_gd_file, out_veg_file = libtts.run_ground_detection(infile = infile, out_gd_file = out_gd_file, out_veg_file = out_veg_file, 
                                                        grid_size=0.1, height_threshold = 0.5)                                                 
"""

import argparse
import numpy as np
from scipy.interpolate import griddata
import cv2

try:
    import matplotlib.pyplot as plt
    PLOTTING_ENABLED = True
except ImportError:
    PLOTTING_ENABLED = False

try:
    from plyfile import PlyData, PlyElement
    PLYFILE_ENABLED = True
except ImportError:
    PLYFILE_ENABLED = False

# --- Core Algorithm Functions ---

def _create_initial_grid(
    points: np.ndarray,
    grid_size: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Internal function to create a raw 2D grid of lowest Z points.
    
    Args:
        points (np.ndarray): Nx3 array of (x, y, z) points.
        grid_size (float): The size of each grid cell.
    
    Returns:
        tuple: A tuple containing:
            - grid_x (np.ndarray): 2D array of X coordinates for the grid.
            - grid_y (np.ndarray): 2D array of Y coordinates for the grid.
            - grid_z_min (np.ndarray): 2D array of minimum Z values for each grid cell.
    """
    points = np.asarray(points)[:, :3]
    min_x, min_y, _ = np.min(points, axis=0)
    max_x, max_y, _ = np.max(points, axis=0)

    x_indices = ((points[:, 0] - min_x) / grid_size).astype(int)
    y_indices = ((points[:, 1] - min_y) / grid_size).astype(int)

    n_x = int((max_x - min_x) / grid_size) + 1
    n_y = int((max_y - min_y) / grid_size) + 1

    grid_z_min = np.full((n_y, n_x), np.inf, dtype=float)
    indices_and_z = np.vstack([y_indices, x_indices, points[:, 2]]).T
    indices_and_z = indices_and_z[indices_and_z[:, 2].argsort()]

    unique_indices, first_occurrence_indices = np.unique(
        indices_and_z[:, :2], axis=0, return_index=True
    )
    unique_indices = unique_indices.astype(int)
    min_z_values = indices_and_z[first_occurrence_indices, 2]

    grid_z_min[unique_indices[:, 0], unique_indices[:, 1]] = min_z_values
    grid_z_min[grid_z_min == np.inf] = np.nan

    grid_x_coords = min_x + np.arange(n_x) * grid_size
    grid_y_coords = min_y + np.arange(n_y) * grid_size
    grid_x, grid_y = np.meshgrid(grid_x_coords, grid_y_coords)

    return grid_x, grid_y, grid_z_min


def _fill_and_smooth_grid(
    grid_z: np.ndarray,
    outlier_std_dev: float,
    gaussian_kernel_size: int
) -> np.ndarray:
    """Internal function to fill, filter, and smooth a 2D grid.

    This helper function performs three main operations on a copy of the input grid:
    1.  Removes statistical outliers by replacing them with NaN.
    2.  Fills any NaN values using linear interpolation based on valid neighbors.
    3.  Applies a Gaussian blur to smooth the entire grid.

    Args:
        grid_z (np.ndarray): A 2D NumPy array representing the Z-values of the grid.
            This grid may contain NaN values.
        outlier_std_dev (float): The number of standard deviations from the mean
            to use as a threshold for outlier removal. Any value outside
            `mean ± outlier_std_dev * std` is replaced with NaN. If this is
            set to 0 or less, this step is skipped.
        gaussian_kernel_size (int): The size of the kernel for Gaussian smoothing.
            This should be a positive, odd integer (e.g., 3, 5, 7). If this
            is set to 0 or less, this step is skipped.

    Returns:
        np.ndarray: A new 2D NumPy array of the same shape as `grid_z` that has
            been filtered, filled, and smoothed.
    """

    processed_grid = np.copy(grid_z)
    
    if outlier_std_dev > 0:
        grid_nonan = processed_grid[~np.isnan(processed_grid)]
        if grid_nonan.size > 0:
            mean_z = np.mean(grid_nonan)
            std_z = np.std(grid_nonan)
            processed_grid[processed_grid > mean_z + outlier_std_dev * std_z] = np.nan
            processed_grid[processed_grid < mean_z - outlier_std_dev * std_z] = np.nan

    nan_mask = np.isnan(processed_grid)
    if np.any(nan_mask):
        valid_coords = np.array(np.nonzero(~nan_mask)).T
        valid_values = processed_grid[~nan_mask]
        
        if valid_coords.size > 0:
            nan_coords = np.array(np.nonzero(nan_mask)).T
            interpolated_values = griddata(valid_coords, valid_values, nan_coords, method='linear')
            # nan_mask_after_linear = np.isnan(interpolated_values)
            # if np.any(nan_mask_after_linear):
            #     nearest_values = griddata(valid_coords, valid_values, nan_coords[nan_mask_after_linear], method='nearest')
            #     interpolated_values[nan_mask_after_linear] = nearest_values
            processed_grid[nan_mask] = interpolated_values

    if gaussian_kernel_size > 0:
        processed_grid = cv2.GaussianBlur(processed_grid, (gaussian_kernel_size, gaussian_kernel_size), 0)

    return processed_grid


def _save_points(
    points_xyzh: np.ndarray,
    outfile: str
):
    """
    A helper function to save point clouds in various formats.

    Args:
        points_xyzh (np.ndarray): An Nx4 array of (x, y, z, h) points.
        outfile (str): The path to the output file. Format is determined by extension.
                       Filename must contain 'xyz', 'xyh', or 'xyzh' to specify columns.
    """
    if outfile is None:
        return

    # Determine which columns to save based on filename
    if "xyzh" in outfile.lower():
        data_to_save = points_xyzh
        dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('h', 'f4')]
    elif "xyh" in outfile.lower():
        data_to_save = points_xyzh[:, [0, 1, 3]]
        dtype = [('x', 'f4'), ('y', 'f4'), ('h', 'f4')]
    else: # Default to xyz
        data_to_save = points_xyzh[:, :3]
        dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4')]

    # Determine file format
    if outfile.lower().endswith(".pts") or outfile.lower().endswith(".txt"):
        np.savetxt(outfile, data_to_save, fmt="%.3f")
        print(f"Saved {len(data_to_save)} points to {outfile}")
    elif outfile.lower().endswith(".ply"):
        if not PLYFILE_ENABLED:
            print("Warning: plyfile library is not installed. Cannot save .ply file.")
            return
        vertex = np.array([tuple(p) for p in data_to_save], dtype=dtype)
        el = PlyElement.describe(vertex, 'vertex')
        PlyData([el]).write(outfile)
        print(f"Saved {len(data_to_save)} points to {outfile}")
    else:
        print(f"Warning: Unknown output file format for {outfile}. Supported formats: .pts, .txt, .ply")

# --- High-Level API Functions ---

def detect_ground(
    points: np.ndarray,
    grid_size: float = 0.5,
    outlier_std_dev: float = 1.0,
    gaussian_kernel_size: int = 5
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Detects the ground surface and returns it as both a point cloud and a grid.

    Args:
        points (np.ndarray): The input Nx3 point cloud.
        grid_size (float): The resolution for the intermediate ground grid.
        outlier_std_dev (float): Std deviations for outlier removal.
        gaussian_kernel_size (int): Kernel size for smoothing.

    Returns:
        A tuple containing:
        - ground_model_points (np.ndarray): Mx3 array of the final ground grid points.
        - grid_x (np.ndarray): 2D array of X coordinates for the grid.
        - grid_y (np.ndarray): 2D array of Y coordinates for the grid.
        - grid_z (np.ndarray): 2D array of Z values for the grid.
    """
    print(f"Creating ground grid with size {grid_size}...")
    grid_x, grid_y, grid_z_initial = _create_initial_grid(points, grid_size)
    
    print("Filling and smoothing ground grid...")
    grid_z_final = _fill_and_smooth_grid(grid_z_initial, outlier_std_dev, gaussian_kernel_size)
    
    # Prepare ground model point cloud for output
    ground_model_points = np.vstack([grid_x.ravel(), grid_y.ravel(), grid_z_final.ravel()]).T
    ground_model_points = ground_model_points[~np.isnan(ground_model_points).any(axis=1)]
    
    return ground_model_points, grid_x, grid_y, grid_z_final


def calculate_height_above_ground(
    points: np.ndarray,
    ground_grid_x: np.ndarray,
    ground_grid_y: np.ndarray,
    ground_grid_z: np.ndarray
) -> np.ndarray:
    """
    Calculates the vertical distance of each point from a pre-computed ground grid.

    Args:
        points (np.ndarray): Nx3 array of (x, y, z) points to normalize.
        ground_grid_x (np.ndarray): 2D array of X coordinates for the ground model.
        ground_grid_y (np.ndarray): 2D array of Y coordinates for the ground model.
        ground_grid_z (np.ndarray): 2D array of Z values for the ground model.

    Returns:
        np.ndarray: An Nx4 array of the original points with an added height column (x, y, z, h).
    """
    if ground_grid_z.size == 0:
        print("Warning: Ground grid is empty. Returning heights as original Z values.")
        return np.column_stack((points, points[:, 2]))

    # Create the set of known points from the grid for interpolation
    grid_points = np.vstack([ground_grid_x.ravel(), ground_grid_y.ravel()]).T
    grid_values = ground_grid_z.ravel()
    
    # Remove any NaN values from the grid before creating the interpolator
    valid_mask = ~np.isnan(grid_values)
    
    if not np.any(valid_mask):
        # If the entire grid is NaN, return zero height for all points
        return np.column_stack((points, np.zeros(points.shape[0])))
        
    # Interpolate the Z value of the ground at the XY location of each point
    interpolated_z = griddata(grid_points[valid_mask], grid_values[valid_mask], points[:, :2], method='linear')
    
    # Fallback for points outside the convex hull of the ground grid
    nan_heights_mask = np.isnan(interpolated_z)
    if np.any(nan_heights_mask):
        nearest_interpolated_z = griddata(grid_points[valid_mask], grid_values[valid_mask], points[nan_heights_mask, :2], method='nearest')
        interpolated_z[nan_heights_mask] = nearest_interpolated_z
        
    heights = points[:, 2] - interpolated_z
    return np.column_stack((points, heights))


def classify_ground_and_vegetation(
    points: np.ndarray,
    height_threshold: float = 0.5,
    out_gd_file: str = None,
    out_veg_file: str = None,
    **ground_detection_params
) -> tuple[np.ndarray, np.ndarray]:
    """
    Performs the full workflow: detect ground, normalize heights, and classify.

    Args:
        points (np.ndarray): The input Nx3 point cloud.
        height_threshold (float): Height to separate vegetation from ground.
        out_gd_file (str, optional): Path to save the ground points.
        out_veg_file (str, optional): Path to save the vegetation points.
        **ground_detection_params: Keyword arguments passed to detect_ground(),
                                   e.g., grid_size=0.5.

    Returns:
        A tuple containing:
        - ground_points (np.ndarray): Nx3 array of points classified as ground.
        - veg_points (np.ndarray): Mx3 array of points classified as vegetation.
    """
    _, grid_x, grid_y, grid_z = detect_ground(points, **ground_detection_params)
    
    normalized_points = calculate_height_above_ground(points, grid_x, grid_y, grid_z)
    
    veg_mask = normalized_points[:, 3] > height_threshold
    ground_mask = ~veg_mask
    
    print(f"Found {np.sum(ground_mask)} ground points and {np.sum(veg_mask)} vegetation points.")
    
    # Save files if paths are provided
    _save_points(normalized_points[ground_mask], out_gd_file)
    _save_points(normalized_points[veg_mask], out_veg_file)
    
    return points[ground_mask], points[veg_mask]

# --- Visualization Functions ---

def plot_ground_model(grid_x: np.ndarray, grid_y: np.ndarray, grid_z: np.ndarray, title: str = "Ground Model (2D Heatmap)"):
    """
    Creates a 2D heatmap directly from the generated ground grid.

    Args:
        grid_x (np.ndarray): 2D array of X coordinates from detect_ground.
        grid_y (np.ndarray): 2D array of Y coordinates from detect_ground.
        grid_z (np.ndarray): 2D array of Z values from detect_ground.
        title (str): The title for the plot.
    """
    if not PLOTTING_ENABLED:
        print("Warning: Matplotlib is not installed. Skipping plot.")
        return
    if grid_z.size == 0:
        print("Warning: Not enough ground model data to create a plot.")
        return

    # Plotting the 2D heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    # Get the extent from the grid coordinates
    extent = [grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()]
    
    im = ax.imshow(grid_z, extent=extent,
                   origin='lower', cmap='terrain') # , interpolation='bilinear'
    
    fig.colorbar(im, ax=ax, label='Z Coordinate (Height)')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_title(title)
    ax.set_aspect('equal', adjustable='box')
    
    print("Displaying 2D ground model plot...")
    plt.show()
    plt.close()

# --- Complete Workflow ---

def run_ground_detection(
    infile: str,
    out_gd_file: str = None,
    out_veg_file: str = None,
    grid_size: float = 0.5,
    height_threshold: float = 0.5,
    **ground_detection_params
) -> list[str,str]: 
    """
    Runs the full ground detection workflow and saves the results.

    Args:
        infile (str): Path to the input point cloud file (.pts or .txt).
        out_gd_file (str, optional): Path to save ground points.
        out_veg_file (str, optional): Path to save vegetation points.
        grid_size (float): Size of the grid cells for DTM generation.
        height_threshold (float): Height threshold to classify vegetation.
        ground_detection_params include:
            - outlier_std_dev (float): Std deviations for outlier removal.
            - gaussian_kernel_size (int): Kernel size for smoothing.

    Returns:
        A tuple containing:
        - Path to the saved ground points file.
        - Path to the saved vegetation points file.
    """
 
    print(f"Loading points from {infile}...")
    if ".pts" in infile.lower() or ".txt" in infile.lower():
        # Load points from .pts or .txt file
        points = np.loadtxt(infile)
    elif ".ply" in infile.lower():
        # Load points from .ply file
        if not PLYFILE_ENABLED:
            raise ImportError("plyfile library is not installed. Cannot load .ply files.")
        ply_data = PlyData.read(infile)
        points = np.array([list(vertex) for vertex in ply_data['vertex'].data])
        # we only need the first three columns (x, y, z)
        points = points[:, :3]
    else:
        raise ValueError("Unsupported file format. Please use .pts, .txt, or .ply files.")

    # Run the full workflow
    classify_ground_and_vegetation(
        points,
        height_threshold=height_threshold,
        out_gd_file=out_gd_file,
        out_veg_file=out_veg_file,
        grid_size=grid_size,
        **ground_detection_params
    )
    
    return out_gd_file, out_veg_file

# --- Main CLI ---

def main():
    """Main function to run ground detection from the command line."""
    parser = argparse.ArgumentParser(description="Detects and classifies ground points from a point cloud.")
    parser.add_argument("infile", help="Path to the input point cloud file (.pts or .txt).")
    parser.add_argument("--out_gd", help="Optional: Path to save ground points.")
    parser.add_argument("--out_veg", help="Optional: Path to save vegetation points.")
    parser.add_argument("--grid_size", type=float, default=0.5, help="Size of the grid cells for DTM generation.")
    parser.add_argument("--height_threshold", type=float, default=0.5, help="Height threshold to classify vegetation.")
    parser.add_argument("--plot", action="store_true", help="Display a 2D plot of the generated ground model.")
    
    args = parser.parse_args()

    try:
        print(f"Loading points from {args.infile}...")
        points = np.loadtxt(args.infile)
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    # --- Run the full workflow ---
    ground_points, veg_points = classify_ground_and_vegetation(
        points,
        height_threshold=args.height_threshold,
        out_gd_file=args.out_gd,
        out_veg_file=args.out_veg,
        grid_size=args.grid_size
    )
    
    # --- Optional Plotting ---
    if args.plot:
        # We need to re-run detect_ground to get the grid for plotting
        _, grid_x, grid_y, grid_z = detect_ground(points, grid_size=args.grid_size)
        plot_ground_model(grid_x, grid_y, grid_z)

    print("Processing complete.")

if __name__ == "__main__":
    main()
