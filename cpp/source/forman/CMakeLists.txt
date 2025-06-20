cmake_minimum_required(VERSION 3.9)
message(STATUS "Forman starts...") #

# xx: add openmp
find_package(OpenMP REQUIRED)

file(GLOB SOURCE_LIBFORMAN "*.h" "*.cpp")
add_library(libFORMAN STATIC ${SOURCE_LIBFORMAN})

TARGET_LINK_LIBRARIES(libFORMAN OpenMP::OpenMP_CXX)

TARGET_LINK_LIBRARIES(libFORMAN libIASTAR)
message(STATUS "Forman ends...") #
