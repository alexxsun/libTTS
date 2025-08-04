## libTTS: A Library for Forest Point Cloud Processing using Topological Data Analysis

**Project Status**
This project is currently under active development. 
The API and features are likely to change frequently, so please check back often for updates.

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

Usage:
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