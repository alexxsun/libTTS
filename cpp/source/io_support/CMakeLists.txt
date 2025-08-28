cmake_minimum_required(VERSION 3.9)
message(STATUS "IOSupport starts...")

# --- Define Search Paths based on Environment ---
set(HEADER_SEARCH_PATHS)
set(LIB_SEARCH_PATHS)

# Check if we are in a Conda environment by looking for its environment variable
if(DEFINED ENV{CONDA_PREFIX})
    message(STATUS "Conda environment detected at: $ENV{CONDA_PREFIX}")
    # Add Conda paths to the FRONT of the search list to give them priority
    list(APPEND HEADER_SEARCH_PATHS "$ENV{CONDA_PREFIX}/include")
    list(APPEND LIB_SEARCH_PATHS "$ENV{CONDA_PREFIX}/lib/LASlib" "$ENV{CONDA_PREFIX}/lib")
else()
    message(STATUS "No Conda environment detected. Searching standard system paths.")
endif()

# Append the standard system paths as fallbacks
list(APPEND HEADER_SEARCH_PATHS "/usr/local/include" "/usr/include")
list(APPEND LIB_SEARCH_PATHS "/usr/local/lib/LASlib" "/usr/local/lib" "/usr/lib")

# show HEADER_SEARCH_PATHS, HEADER_SEARCH_PATHS
message(STATUS "HEADER_SEARCH_PATHS: ${HEADER_SEARCH_PATHS}")
message(STATUS "LIB_SEARCH_PATHS: ${LIB_SEARCH_PATHS}")

# --- Find LASlib ---
# Find the path to the main header file using the prioritized search paths
find_path(LASLIB_INCLUDE_DIR
    NAMES LASlib/lasreader.hpp
    HINTS ${HEADER_SEARCH_PATHS}
    DOC "Path to the LASlib header directory"
)

# Find the compiled library file using the prioritized search paths
find_library(LASLIB_LIBRARY
    NAMES LASlib
    HINTS ${LIB_SEARCH_PATHS}
    DOC "Path to the LASlib library file"
)

# Check if both were found
if(NOT LASLIB_INCLUDE_DIR OR NOT LASLIB_LIBRARY)
    message(FATAL_ERROR "Could not find the installed LASlib library. Make sure it is installed in /usr/local/ or a standard system path.")
else()
    message(STATUS "Found LASlib headers: ${LASLIB_INCLUDE_DIR}")
    message(STATUS "Found LASlib library: ${LASLIB_LIBRARY}")
endif()

# --- Create the Interface Library ---
add_library(libIOSUPPORT INTERFACE)

# Add the local directory for happly.h
target_include_directories(libIOSUPPORT
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add the found LASlib include directory
target_include_directories(libIOSUPPORT
    INTERFACE
        ${LASLIB_INCLUDE_DIR}
)

# Link against the found LASlib library
target_link_libraries(libIOSUPPORT
    INTERFACE
        ${LASLIB_LIBRARY}
)

message(STATUS "IOSupport ends...")
