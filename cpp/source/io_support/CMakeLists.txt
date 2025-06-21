cmake_minimum_required(VERSION 3.9)
message(STATUS "IOSupport starts...")

file(GLOB SOURCE_IOSUPPORT "*.h" "*.cpp")
add_library(libIOSUPPORT STATIC ${SOURCE_IOSUPPORT})
message(STATUS "IOSupport ends...")

