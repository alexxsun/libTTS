# python/src/libtts/__init__.py

# When scikit-build-core compiles the C++ code, it creates a module file
# (e.g., _libtts.cpython-311-x86_64-linux-gnu.so) 
# check cpp/cmakelists.txt and python/pyproject.toml to define the module name and install path.
# (e.g., site-packages/libtts/_libtts.cpython-311-x86_64-linux-gnu.so)
# This module is imported here to make its functions available in the package namespace.
from ._libtts import (generate_alpha_shape_cpp as generate_alpha_shape,
                      get_oversegments_cpp as get_oversegments,
                      tls_extract_single_trees_cpp as tls_extract_single_trees,
                      als_segment as als_extract_single_trees)

# Import any high-level Python functions or classes to expose.
from .ground_detection import detect_ground, calculate_height_above_ground, classify_ground_and_vegetation, plot_ground_model, run_ground_detection
from .tree_detection import filter_tree_bases, cluster_points_dbscan, run_tree_detection
from .points_downsampling import (downsample_by_lastools, downsample_object_based, 
                                  downsample_object_based_from_mesh, downsample_object_based_from_points, 
                                  run_downsampling)
from .label_propagation import label_points_layered_nn, label_points_region_growing, run_label_propagation

# Define the package version in the top-level __init__.py.
# This allows users to check the version with `libtts.__version__`.
__version__ = "0.0.1"

# The __all__ variable defines the public API in the package. 
# It's a listof strings containing the names of the objects to be imported 
# when a user does `from libtts import *`. 
# It is good practice to define this.
__all__ = [
    # Exposing the C++ core functions
    "generate_alpha_shape",
    "get_oversegments",
    "tls_extract_single_trees",
    "als_extract_single_trees",
    # Exposing the Python functions
    "detect_ground",
    "calculate_height_above_ground",
    "classify_ground_and_vegetation",
    "plot_ground_model",
    "run_ground_detection",
    "filter_tree_bases",
    "cluster_points_dbscan",
    "run_tree_detection",
    ### Points downsampling
    "downsample_by_lastools",
    "downsample_object_based",
    "downsample_object_based_from_mesh",
    "downsample_object_based_from_points",
    "run_downsampling",
    ### Label propagation
    "label_points_layered_nn",
    "label_points_region_growing",
    "run_label_propagation",
    # "common",
    "__version__"
]
