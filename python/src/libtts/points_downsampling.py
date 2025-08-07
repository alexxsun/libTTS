"""Functions for downsampling point clouds using two distinct methods.

This module provides functionality to reduce the number of points in a point
cloud file. It supports two primary strategies:

1.  **LAStools-based**: A general-purpose random downsampling method that
    uses the external LAStools 'las2las' executable.
2.  **Object-based**: An algorithm designed for wood-leaf separation in
    segmented point clouds. It filters points within each segment based on
    local point density, effectively removing sparse "leaf" points while
    retaining dense "woody" structures.

Dependencies:
    - numpy
    - plyfile
    - laspy (for LAStools-based method)
    - scikit-learn (for object-based method)

External Dependencies:
    - LAStools: Required for the 'lastools' method. The path to the 'bin'
      directory must be provided.

Example:
    .. code-block:: python

        # LAStools-based downsampling to keep 50% of points
        ds_file = libtts.run_downsampling(
            method="lastools",
            infile="input.ply",
            lastools_bin="/path/to/lastools/bin/",
            fraction=0.5
        )

        # Object-based downsampling on a pre-segmented point cloud
        woody_pts_file = libtts.run_downsampling(
            method="object_based",
            infile="oversegmented_tree.pts",
            input_type="overseg",
            distance_threshold=0.1
        )

        # Object-based downsampling starting from a raw point cloud (requires C++ module)
        woody_pts_file_from_raw = libtts.run_downsampling(
            method="object_based",
            infile="raw_tree.pts",
            input_type="pts",
            alpha_sq=0.01,
            distance_threshold=0.1
        )
"""

import argparse
import multiprocessing as mp
import numpy as np
import subprocess
import os
import sys
import tempfile

# --- Optional Imports ---
# These are encapsulated in functions or checked to provide clear error messages.
try:
    from plyfile import PlyData, PlyElement
except ImportError:
    print("Error: The 'plyfile' library is required. Please run 'pip install plyfile'")
    sys.exit(1)

try:
    import laspy
except ImportError:
    # This is only a problem if using the 'lastools' method.
    # A check is performed in the relevant function.
    pass

try:
    from sklearn.neighbors import NearestNeighbors
except ImportError:
    # This is only a problem if using the 'object_based' method.
    # A check is performed in the relevant function.
    pass

# --- Optional C++ Module for Object-Based Method ---
CPP_MODULE_AVAILABLE = True
try:
    from ._libtts import get_oversegments_cpp as _oversegment_tree_cpp
    from ._libtts import generate_alpha_shape_cpp as _alpha_shape_cpp
except ImportError:
    print("Info: C++ module 'libtts' not found. Advanced object-based features will be unavailable.")
    CPP_MODULE_AVAILABLE = False


# #############################################################################
#
# LAStools BASED DOWNSAMPLING
#
# #############################################################################

def _ply_to_las(ply_path: str, las_path: str):
    """Converts a PLY file's XYZ data to a LAS 1.2 file.

    Args:
        ply_path (str): Path to the input .PLY file.
        las_path (str): Path for the output .LAS file.
    """
    print(f"Converting {ply_path} to LAS...")
    plydata = PlyData.read(ply_path)
    vertices = plydata['vertex']
    points = np.vstack([vertices['x'], vertices['y'], vertices['z']]).T

    header = laspy.LasHeader(version="1.2", point_format=0)
    header.offsets = np.min(points, axis=0)
    header.scales = np.array([0.001, 0.001, 0.001])

    las = laspy.LasData(header)
    las.x = points[:, 0]
    las.y = points[:, 1]
    las.z = points[:, 2]

    las.write(las_path)
    print(f"Successfully created temporary LAS file: {las_path}")


def _las_to_ply(las_path: str, ply_path: str):
    """Converts a LAS file's XYZ data to a PLY file.

    Args:
        las_path (str): Path to the input .LAS file.
        ply_path (str): Path for the output .PLY file.
    """
    print(f"Converting {las_path} back to PLY...")
    las = laspy.read(las_path)
    points = np.vstack((las.x, las.y, las.z)).transpose()

    vertex = np.array([tuple(p) for p in points], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4')])
    el = PlyElement.describe(vertex, 'vertex')

    PlyData([el], text=False).write(ply_path)
    print(f"Successfully created final PLY file: {ply_path}")


def downsample_by_lastools(infile: str, lastools_bin_dir: str, keep_random_fraction: float) -> str:
    """Downsamples a point cloud file by a random fraction using LAStools.

    This method performs a multi-step process:

    1. Converts the input file (.ply, .pts) to a temporary LAS file if needed.
    2. Calls the `las2las` executable to perform random downsampling.
    3. Converts the downsampled LAS file back to a PLY file for the final output.

    Args:
        infile (str): Path to the input file (.ply, .las, or .pts).
        lastools_bin_dir (str): Path to the LAStools 'bin' directory.
        keep_random_fraction (float): The fraction of points to keep (e.g., 0.1 for 10%).

    Returns:
        str: The path to the final, downsampled .PLY output file.

    Raises:
        ImportError: If the 'laspy' library is not installed.
        FileNotFoundError: If the 'las2las' or 'txt2las' executable is not found.
        RuntimeError: If a LAStools command fails.
        ValueError: If the input file format is unsupported.
        NotImplementedError: If run on a non-Linux environment.
    """
    if 'laspy' not in sys.modules:
        raise ImportError("The 'laspy' library is required for the LAStools method. Please run 'pip install laspy'")

    las2las_exe = os.path.join(lastools_bin_dir, "las2las")
    if sys.platform == "win32":
        # we don't support windows yet
        raise NotImplementedError("Please use a Linux environment.")

    if not os.path.exists(las2las_exe):
        raise FileNotFoundError(f"las2las executable not found at: {las2las_exe}")

    base_name = os.path.splitext(os.path.basename(infile))[0]
    output_dir = os.path.dirname(infile)
    outfile = os.path.join(output_dir, f"{base_name}_ds_{keep_random_fraction}.ply")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_las_in = os.path.join(temp_dir, "temp_in.las")
        temp_las_out = os.path.join(temp_dir, "temp_out.las")

        if infile.endswith('.ply'):
            _ply_to_las(infile, temp_las_in)
        elif infile.endswith('.las'):
            # If the input is already a LAS file, just copy it to the temp directory
            temp_las_in = infile
        elif infile.endswith('.pts'):
            # convert PTS to LAS
            txt2las_exe = os.path.join(lastools_bin_dir, "txt2las")
            command = [
                txt2las_exe, "-i", infile, "-o", temp_las_in, "-parse", "xyz"
            ]
            try:
                subprocess.run(command, check=True, capture_output=True, text=True, shell=False)
            except subprocess.CalledProcessError as e:
                print("--- LAStools Error ---")
                print(e.stderr)
                print("----------------------")
                raise RuntimeError("las2las command failed for PTS to LAS conversion.") from e
        else:
            raise ValueError(f"Unsupported file format: {infile}. Supported formats are .ply, .las, and .pts.")

        print(f"Running las2las to keep {keep_random_fraction * 100:.1f}% of points...")
        command = [las2las_exe, "-i", temp_las_in, "-o", temp_las_out, "-keep_random_fraction", str(keep_random_fraction)]

        try:
            subprocess.run(command, check=True, capture_output=True, text=True, shell=False)
        except subprocess.CalledProcessError as e:
            print("--- LAStools Error ---")
            print(e.stderr)
            print("----------------------")
            raise RuntimeError("las2las command failed.") from e

        _las_to_ply(temp_las_out, outfile)

    return outfile


# #############################################################################
#
# OBJECT-BASED DOWNSAMPLING (WOOD-LEAF SEPARATION)
#
# #############################################################################

def _select_pts_per_label(points_in_segment: np.ndarray, th_neighbors: int, th_distance: float) -> np.ndarray:
    """Filters points in a segment based on average neighbor distance.

    For a given array of points (assumed to be from the same segment), this
    function calculates the k-nearest neighbors for each point. It keeps
    only those points whose average distance to its neighbors is below a
    specified threshold, effectively filtering for dense clusters.

    Args:
        points_in_segment (np.ndarray): An array of points (NxM, M>=4) belonging to one segment.
        th_neighbors (int): The number of neighbors to consider for the distance calculation.
        th_distance (float): The mean distance threshold. Points with a mean neighbor
                             distance below this value will be kept.

    Returns:
        np.ndarray: A NumPy array of points that met the distance criteria. Returns
        an empty array if the input segment has fewer points than `th_neighbors`.
    """
    if len(points_in_segment) < th_neighbors:
        return np.array([])

    xyz_points = points_in_segment[:, :3]
    nbrs = NearestNeighbors(n_neighbors=th_neighbors, algorithm='ball_tree').fit(xyz_points)
    distances, _ = nbrs.kneighbors(xyz_points)

    mean_distances = np.mean(distances[:, 1:], axis=1)

    selected_indices = np.where(mean_distances < th_distance)[0]
    return points_in_segment[selected_indices]


def downsample_object_based(infile: str, th_avg_dis: float, num_processes: int = 8) -> str:
    """Filters a labeled point cloud to separate woody from leafy points.

    This function reads a point cloud with labels (e.g., from an oversegmentation
    result), groups points by their label, and then filters each group in
    parallel. The filtering keeps points in dense clusters, effectively
    removing sparse "leaf" points.

    Args:
        infile (str): Path to the input point cloud file (.pts or .ply with labels).
                      The file must contain XYZ and label data.
        th_avg_dis (float): The average distance threshold for classifying a point
                            as "woody".
        num_processes (int): The number of CPU cores to use for parallel processing.

    Returns:
        str: The path to the output file containing the filtered "woody" points.

    Raises:
        ImportError: If 'scikit-learn' is not installed.
        ValueError: If the input file format is not supported, if it lacks a 'label'
                    property, or if no points are selected after filtering.
        FileNotFoundError: If the input file does not exist.
    """
    if 'sklearn' not in sys.modules:
        raise ImportError("The 'scikit-learn' library is required for the object-based method. Please run 'pip install scikit-learn'")

    if not os.path.exists(infile):
        raise FileNotFoundError(f"Input file not found: {infile}")

    file_ext = os.path.splitext(infile)[1]
    if file_ext not in ['.pts', '.ply']:
        raise ValueError(f"File extension '{file_ext}' is not supported for object-based method. Use '.pts' or '.ply'.")

    print("Reading input file...")
    if file_ext == '.pts':
        pts_xyzl = np.loadtxt(infile)
    elif file_ext == '.ply':
        plydata = PlyData.read(infile)
        if 'label' not in plydata['vertex'].data.dtype.names:
            raise ValueError("Input .ply file for object-based method must contain a 'label' property.")
        pts_xyzl = np.array([
            [x, y, z, l] for x, y, z, l in zip(
                plydata['vertex']['x'], plydata['vertex']['y'], plydata['vertex']['z'], plydata['vertex']['label']
            )
        ])

    print(f"Total points read: {len(pts_xyzl)}")

    grouped_points = {}
    for point in pts_xyzl:
        label = point[3]
        if label not in grouped_points:
            grouped_points[label] = []
        grouped_points[label].append(point)

    segments = [np.array(points) for points in grouped_points.values()]
    print(f"Number of unique segments: {len(segments)}")

    print(f"Starting parallel processing on {num_processes} cores...")
    pool_args = [(segment, 4, th_avg_dis) for segment in segments]
    with mp.Pool(processes=num_processes) as pool:
        results = pool.starmap(_select_pts_per_label, pool_args)

    selected_pts_list = [res for res in results if res.size > 0]
    if not selected_pts_list:
        raise ValueError("No points were selected after filtering. Try increasing the distance threshold.")

    selected_pts = np.vstack(selected_pts_list)
    print(f"Selected 'woody' points: {len(selected_pts)}")

    outfile = f"{os.path.splitext(infile)[0]}_woody{file_ext}"
    print(f"Saving output to: {outfile}")

    if file_ext == ".pts":
        np.savetxt(outfile, selected_pts, fmt="%.3f")
    elif file_ext == ".ply":
        # Create a structured numpy array for plyfile
        vertex_data = np.array(
            [tuple(row) for row in selected_pts],
            dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('label', 'i4')]
        )
        el = PlyElement.describe(vertex_data, 'vertex')
        PlyData([el], text=False).write(outfile)

    return outfile

def downsample_object_based_from_mesh(infile, th_avg_dis=0.1):
    """Runs oversegmentation on a mesh file, then downsamples the result.

    This function provides a workflow for mesh inputs. It first calls a C++
    backend function to perform oversegmentation on the mesh, then passes
    the resulting labeled point cloud to the object-based downsampler.

    Args:
        infile (str): Path to the input mesh file.
        th_avg_dis (float): The average distance threshold passed to the
                            final downsampling step. Defaults to 0.1.

    Returns:
        str: The path to the final "woody" points file.

    Raises:
        NotImplementedError: If the required C++ module is not available.
    """
    # use cpp function to process infile
    if not CPP_MODULE_AVAILABLE:
        raise NotImplementedError(
            "The C++ module 'libtts' is required for this method. "
            "Please ensure it is compiled and available in your Python environment."
        )

    # Call the C++ function
    overseg_file = _oversegment_tree_cpp(infile) # infile: mesh file

    # process the overseg_file
    print(f"Oversegmentation results saved to {overseg_file}")
    return downsample_object_based(overseg_file, th_avg_dis)

def downsample_object_based_from_points(infile, th_alpha_sq = 0.01, th_avg_dis=0.1):
    """Runs a full wood-leaf separation workflow from a raw point cloud.

    This function orchestrates a multi-step C++ pipeline:

    1. Generates an alpha shape from the raw points.
    2. Performs oversegmentation on the resulting mesh.
    3. Passes the labeled point cloud to the object-based downsampler.

    Args:
        infile (str): Path to the input raw point cloud file (.pts).
        th_alpha_sq (float): Alpha squared value for alpha shape generation.
                             Defaults to 0.01.
        th_avg_dis (float): The average distance threshold for the final
                            downsampling step. Defaults to 0.1.

    Returns:
        str: The path to the final "woody" points file.

    Raises:
        NotImplementedError: If the required C++ module is not available.
    """
    # use cpp function to process infile
    if not CPP_MODULE_AVAILABLE:
        raise NotImplementedError(
            "The C++ module 'libtts' is required for this method. "
            "Please ensure it is compiled and available in your Python environment."
        )

    # Call the C++ function
    as_file = _alpha_shape_cpp(infile, th_alpha_sq)
    overseg_file = _oversegment_tree_cpp(as_file) # infile: mesh file

    # process the overseg_file
    print(f"Oversegmentation results saved to {overseg_file}")
    return downsample_object_based(overseg_file, th_avg_dis)

def run_downsampling(method: str, **kwargs):
    """High-level wrapper to run the selected downsampling method.

    This function acts as a dispatcher, calling the appropriate low-level
    downsampling function based on the specified method and its arguments.

    Args:
        method (str): The downsampling method to use ('lastools' or 'object_based').
        **kwargs: Keyword arguments for the underlying methods.
            For 'lastools', requires: `infile` (str), `lastools_bin` (str),
            and `fraction` (float).
            For 'object_based', requires: `infile` (str), `input_type` (str),
            and `th_avg_dis` (float).

    Returns:
        str: The path to the output file.

    Raises:
        ValueError: If the method is unknown or required arguments are missing.
    """
    infile = kwargs.get("infile")
    if not infile:
        raise ValueError("An input file ('infile') must be provided.")

    if method == "lastools":
        if 'lastools_bin' not in kwargs or 'fraction' not in kwargs:
            raise ValueError("For 'lastools' method, 'lastools_bin' and 'fraction' are required.")
        return downsample_by_lastools(infile, kwargs['lastools_bin'], kwargs['fraction'])

    elif method == "object_based":
        input_type = kwargs.get('input_type', 'overseg')
        th_avg_dis = kwargs.get('distance_threshold', 0.1)
        num_processes = kwargs.get('processes', 8)

        if input_type == "overseg":
            # This directly calls the function that handles labeled point clouds
            return downsample_object_based(infile, th_avg_dis, num_processes)
        elif input_type == "mesh":
            # This calls the function that processes mesh files
            return downsample_object_based_from_mesh(infile, th_avg_dis)
        elif input_type == "pts":
            # This calls the function that processes point cloud files
            return downsample_object_based_from_points(infile, kwargs.get('alpha_sq', 0.01), th_avg_dis)
        else:
            raise ValueError(f"Unknown input_type for object_based method: {input_type}")
    else:
        raise ValueError(f"Unknown method: {method}")


# #############################################################################
#
# MAIN COMMAND-LINE INTERFACE
#
# #############################################################################

def main():
    """Main function to handle command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Downsample a point cloud file using LAStools or an object-based method.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="method", required=True, help="The downsampling method to use.")

    # --- LAStools Subparser ---
    parser_lastools = subparsers.add_parser(
        "lastools",
        help="Downsample by keeping a random fraction of points using LAStools.",
        description="""
Method: LAStools
Downsamples a .PLY file by converting it to LAS, running the external las2las
executable to keep a random fraction of points, and converting it back to PLY.
Requires a local LAStools installation.
"""
    )
    parser_lastools.add_argument("-i", "--infile", required=True, help="Path to the input .PLY file.")
    parser_lastools.add_argument("-l", "--lastools_bin", required=True, help="Path to the LAStools 'bin' directory.")
    parser_lastools.add_argument("-f", "--fraction", type=float, required=True, help="Fraction of points to keep (e.g., 0.1 for 10%%).")

    # --- Object-Based Subparser ---
    parser_object = subparsers.add_parser(
        "object_based",
        help="Filter a labeled point cloud to separate wood and leaf components.",
        description="""
Method: Object-Based (Wood-Leaf Separation)
Filters a labeled point cloud (.pts or .ply with a 'label' property).
It assumes the input is an oversegmentation result where each segment
represents a part of an object (e.g., a tree). It removes sparse points
(leaves) by keeping only points in dense areas within each segment.
"""
    )
    parser_object.add_argument("-i", "--infile", required=True, help="Path to the input labeled file (.pts or .ply).")
    parser_object.add_argument(
        "-d", "--distance_threshold", type=float, default=0.1,
        help="The average distance threshold for keeping a point. Default is 0.1."
    )
    parser_object.add_argument(
        "-p", "--processes", type=int, default=8,
        help="Number of CPU cores to use for parallel processing. Default is 8."
    )
    parser_object.add_argument(
        "-t", "--input_type", type=str, default="overseg", choices=["overseg", "mesh", "pts"],
        help="Type of input data. 'mesh' and 'pts' require a compiled C++ module. Default: 'overseg'."
    )


    args = parser.parse_args()
    
    try:
        # Use the run_downsampling wrapper function to handle the logic
        kwargs = vars(args)
        method = kwargs.pop('method')
        output_file = run_downsampling(method=method, **kwargs)
        
        if output_file:
            print(f"\nDownsampling complete. Output saved to: {output_file}")
        else:
            print("\nOperation finished, but no output file was generated (e.g., requires C++ module).")

    except (FileNotFoundError, RuntimeError, ValueError, ImportError, NotImplementedError, Exception) as e:
        print(f"\nAn error occurred: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
