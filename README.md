## libTTS: A Library for Forest Point Cloud Processing using Topological Data Analysis

**Project Status**
This project is currently under active development. 
The API and features may change in the near future.

### Cpp library
You can compile the library from source and have the executables `xx_tts` for processing point clouds.

Compilation requires CMake and a C++ compiler (e.g., g++ or clang++).
```
cd libtts_project_folder
mkdir build
cd build
cmake ../
make -j4
```

Usage:
```
xx_tts [options] <input_file> <output_file>
```
### Python package 

Install in editable mode for development:
```
cd libtts_project_folder/python
pip install -e .[dev]
pip install --verbose -e .
```

Uninstall 
```
pip uninstall libtts
```

Usage:
```python
from libtts import xxx
```
More examples are available in the `examples/` directory of the Python package.


### Project Structure
TOC

```
libtts_project/
├── cpp/                     # All C++ code lives here
│   ├── src/                 # Your core C++ library source files (.cpp, .h)
│   ├── bindings/            # C++ code specifically for Python bindings
│   │   └── bindings.cpp
│   └── CMakeLists.txt       # CMake file for building ONLY the C++ part
│
├── python/                  # All Python code lives here
│   ├── pyproject.toml       # Configuration for the Python package build
│   └── src/
│       └── libtts/          # The Python package source
│           ├── __init__.py
│           └── analysis.py
│
├── some_examples/
|
└── README.md
```

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