cmake_minimum_required(VERSION 3.9)
message(STATUS "IOSupport starts...")

# --- Find the LASlib package ---
# This command will find LASlib whether it's in a Conda prefix or a system path.
find_package(LASlib REQUIRED)
message(STATUS "Found LASlib package")

# --- Set the RPATH for the final installed library ---
# Get the full path to the imported LASlib shared library from the package.
get_target_property(LASLIB_LOCATION LASlib IMPORTED_LOCATION_RELEASE)
if(NOT LASLIB_LOCATION)
    # Fallback for non-Release or single-configuration builds
    get_target_property(LASLIB_LOCATION LASlib IMPORTED_LOCATION)
endif()

if(LASLIB_LOCATION)
    # Get the directory containing the library
    get_filename_component(LASLIB_LIB_DIR ${LASLIB_LOCATION} DIRECTORY)
    message(STATUS "Setting install RPATH to LASlib directory: ${LASLIB_LIB_DIR}")

    # Set the RPATH that will be baked into the final installed .so file.
    # This tells the OS where to find libLASlib.so at runtime.
    # We set this as a CACHE variable to ensure it's passed up to the main project.
    set(CMAKE_INSTALL_RPATH "${LASLIB_LIB_DIR}" CACHE STRING "Runtime path for final library" FORCE)
else()
    message(WARNING "Could not determine LASlib location to set RPATH.")
endif()

# --- Create the Interface Library ---
add_library(libIOSUPPORT INTERFACE)

# Add the local directory for happly.h
target_include_directories(libIOSUPPORT
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link against the imported target provided by the LASlib package.
# This single command correctly sets up include directories and the library link.
# As defined in the package's .cmake files, the target is named "LASlib".
target_link_libraries(libIOSUPPORT
    INTERFACE
        LASlib
)

message(STATUS "IOSupport ends...")
