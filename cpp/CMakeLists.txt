# A higher version allows for more modern CMake features.
cmake_minimum_required(VERSION 3.15)

# Define the project name and specify CXX (C++) as the language.
project(libtts_core LANGUAGES CXX)

# --- Compiler and Build Configuration ---
# Set the C++ standard to C++17.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set Position Independent Code (PIC) for all targets.
# This is required for shared libraries and Python bindings.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- Build C++ Libraries ---
add_subdirectory(source)

# --- Build Standalone Executable ---
option(BUILD_EXECUTABLE "Build the standalone C++ executable" ON)
if(BUILD_EXECUTABLE)
    message(STATUS "Building standalone executable 'xx_tts'...")
    add_executable(xx_tts source/tts_test_main.cpp)

    # Link all required libraries.
    target_link_libraries(xx_tts PRIVATE
        libFORMAN
        libProjects
        libUSAGE
        libIOSUPPORT
    )
endif()


# --- Build Python Bindings (Conditionally) ---
# This is the main switch controlled by the pyproject.toml.
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(BUILD_PYTHON_BINDINGS)
    message(STATUS "Building Python bindings...")

    # Find the pybind11 package, which is provided by scikit-build-core.
    find_package(pybind11 REQUIRED)

    # Create the Python extension module named '_libtts'.
    pybind11_add_module(_libtts
        SHARED
        bindings/bindings.cpp
    )

    # Link the Python module against the necessary C++ libraries.
    target_link_libraries(_libtts PRIVATE libProjects)

    message(STATUS "Python bindings configuration complete.")
endif()

## show some values for debugging
message(STATUS ": ${SKBUILD}")
message(STATUS "${Python_SITELIB}")
