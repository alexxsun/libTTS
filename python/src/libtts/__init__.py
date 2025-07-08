# python/src/libtts/__init__.py

# When scikit-build-core compiles the C++ code, it creates a module file
# (e.g., _libtts.cpython-311-x86_64-linux-gnu.so) 
# check cpp/cmakelist.txt and python/pyproject.toml to define the module name and install path.
# (e.g., site-packages/libtts/_libtts.cpython-311-x86_64-linux-gnu.so)
# This module is imported here to make its functions available in the package namespace.
from ._libtts import alpha_shape_generation, get_oversegments, extract_single_trees, als_segment

# Import any high-level Python functions or classes to expose.
from .simple_ground_dection import normalize_pts, classify_pts
from .simple_tree_detection import detect_trees
from .object_downsampling import extract_woody_points, downsample_points
from .label_propagation import label_points_layered_nn, label_points_region_growing#as propagate_labels

# Define the package version in the top-level __init__.py.
# This allows users to check the version with `libtts.__version__`.
__version__ = "0.0.1"

# The __all__ variable defines the public API in the package. 
# It's a listof strings containing the names of the objects you want to be imported 
# when a user does `from libtts import *`. 
# It is good practice to define this.
__all__ = [
    # Exposing the C++ core functions
    "alpha_shape_generation",
    "get_oversegments",
    "extract_single_trees",
    "als_segment",
    # Exposing the Python functionss
    "normalize_pts",
    "classify_pts",
    "detect_trees",
    "extract_woody_points",
    "downsample_points",
    "label_points_layered_nn",
    "label_points_region_growing"
]
