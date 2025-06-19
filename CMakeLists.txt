cmake_minimum_required(VERSION 3.9) # cmake 3.9 has better support for openmp
message(STATUS "Start running cmake...")

project(TTS CXX)
enable_language(CXX)
# C++17 std::filesystem
# todo: use C++ 20, support string format
# Set C++17 standard and enable optimization for Release build
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "Configuring for Release mode")
set(CMAKE_BUILD_TYPE Release)

add_subdirectory(source)

# tree segmentation.
message(STATUS "tts starts...") #
add_executable(xx_tts "source/tts_test_main.cpp")
TARGET_LINK_LIBRARIES(xx_tts libFORMAN)
TARGET_LINK_LIBRARIES(xx_tts libProjects)
TARGET_LINK_LIBRARIES(xx_tts libUSAGE)
message(STATUS "tts done")

# --- pybind11 ---
set (PYTHON_EXECUTABLE "/home/alex/.pyenv/shims/python")
# Position Independent Code (PIC)
set_target_properties(libIASTAR PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(libFORMAN PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(libANN PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(libProjects PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(libUSAGE PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Find pybind11.
find_package(pybind11 REQUIRED)

# --- Python Bindings ---
# pybind11_add_module will create a Python extension module
pybind11_add_module(xx_tts_py             # CMake target name for the module
                    SHARED                # Build as a shared library for Python to import
                    bindings.cpp          # The binding code
)

# Link your C++ library to the bindings
target_link_libraries(xx_tts_py PRIVATE libProjects)

message(STATUS "pybinding done")