cmake_minimum_required(VERSION 3.9)
message(STATUS "IOSupport starts...")

#file(GLOB SOURCE_IOSUPPORT "*.h" "*.cpp")
#add_library(libIOSUPPORT STATIC ${SOURCE_IOSUPPORT})

# if only .h file exist in this folder, it is a header-only library
# Create an INTERFACE library target. This doesn't compile anything.
add_library(libIOSUPPORT INTERFACE)

# Specify that any target linking to libIOSUPPORT should have the
# current directory added to its include path.
target_include_directories(libIOSUPPORT
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
)

message(STATUS "IOSupport ends...")

