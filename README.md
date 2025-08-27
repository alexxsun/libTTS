## libTTS: A Library for Forest Point Cloud Processing using Topological Data Analysis

**Project Status**
This project is currently under active development. 
The API and features are likely to change frequently, so please check back often for updates.

### Quick Start

1. Set up Conda Environment
```
# Create a new conda environment (e.g., named 'libtts_env') with Python 3.12
$ conda create -n libtts_env python=3.12
# or if you want to create it in a specific path
$ conda create -p /path/to/your/pyenvs/libtts_env python=3.12
# Activate the new environment
$ conda activate libtts_env
# Install CGAL and other C++ dependencies from conda-forge
$ conda install -c conda-forge cmake boost cgal-cpp cxx-compiler
```

2. Clone the Repository
```
$ git clone https://github.com/alexxsun/libTTS
$ cd libTTS
```

3. Install the Python Package
```
$ cd python/
$ pip install .
```

4. Run a Quick Test
```
$ cd ../some_examples/
$ python quick_test.py
```

5. Clean Up (Optional)
```
# Deactivate the current environment
$ conda deactivate
# Remove the environment and all its packages
$ conda env remove -n libtts_env
```


### Dependencies

C++:
- C++ compiler (e.g., g++, clang++)
- CMake 
- Boost 
- CGAL 
- openmp
- LASlib from LAStools

Python:
- check the `pyproject.toml` file for the Python package dependencies.

Bindings:
- pybind11 (for C++ to Python bindings)


### Cpp library

CGAL library is required for the C++ library.

You can compile the library from source and have the executables `xx_tts` for processing point clouds.
Compilation requires CMake and a C++ compiler (e.g., g++ or clang++). 
```
cd libtts_project_folder/cpp
mkdir build
cd build
cmake ../
make -j4
```

Example usage:
```
xx_tts infile.ply/.pts alpha_value_sq -as2
xx_tts inmesh.ply/.off tree_locations.pts -tts
```

### Python package 

Install in editable mode for development:
```
cd libtts_project_folder/python
pip install --verbose -e .
```

Uninstall 
```
pip uninstall libtts
```

Usage:
Please check examples in the `some_examples/` folder.


### Project Structure
Source code

```
libtts_project/
├── cpp/                     # Top-level directory for all C++ code
│   ├── bindings/            # C++ code for Python bindings (pybind11)
│   │   └── bindings.cpp
│   ├── CMakeLists.txt       # Build script for the C++ library
│   ├── data/                # Needed C++ data files
│   └── source/              # Core C++ library source files (.cpp, .h)
├── python/                  # Top-level directory for the Python package
│   ├── docs/                # Sphinx documentation source
│   │   ├── build/           # Compiled documentation (e.g., HTML files)
│   │   └── source/          # Sphinx configuration and .rst files
│   │       ├── api.rst
│   │       ├── conf.py
│   │       └── index.rst
│   ├── pyproject.toml       # Python package build configuration
│   ├── src/                 # Main Python source code directory
│   │   └── libtts/          # The Python package itself
│   │       ├── __init__.py
│   │       ├── ground_detection.py
│   │       ├── label_propagation.py
│   │       ├── points_downsampling.py
│   │       └── tree_detection.py
│   └── tests/               #
├── some_examples/           # Example scripts and data
├── README.md                # This file
└── LICENSE                  # Project license file
```

After the installation using `pip`,  the package and files will be located at in the Python `site-packages` folder.

### Ducumentation
Local documentation is available in the `docs/` folder.

Need to install Sphinx and the necessary dependencies:
- sphinx
- sphinx-rtd-theme

And run
```
cd libtts_project_folder/docs
make html
```
to generate the HTML documentation.
The generated documentation will be available in `docs/build/html/index.html`.

You can also view the documentation online at [libTTS Documentation](https://libtts.readthedocs.io/en/latest/).


### Citation
The methods implemented in this library are described in the following publications. 

For Terrestrial Laser Scanning (TLS) Data:
```
@article{xu2023topology_tls,
  title={Topology-based individual tree segmentation for automated processing of terrestrial laser scanning point clouds},
  author={Xu, Xin and Iuricich, Federico and Calders, Kim and Armston, John and De Floriani, Leila},
  journal={International Journal of Applied Earth Observation and Geoinformation},
  volume={116},
  pages={103145},
  year={2023},
  publisher={Elsevier}
}
```

For Airborne Laser Scanning (ALS) Data:
```
@article{xu2023topology_als,
  title={A topology-based approach to individual tree segmentation from airborne LiDAR data},
  author={Xu, Xin and Iuricich, Federico and De Floriani, Leila},
  journal={GeoInformatica},
  pages={1--30},
  year={2023},
  publisher={Springer}
}
```