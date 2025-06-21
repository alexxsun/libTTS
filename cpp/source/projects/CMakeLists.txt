cmake_minimum_required(VERSION 3.9)
message(STATUS "Projects starts...") #
file(GLOB SOURCE_LIBProjects "*.h" "*.cpp")
add_library(libProjects STATIC ${SOURCE_LIBProjects})

find_package(CGAL REQUIRED)

target_link_libraries(libProjects CGAL::CGAL)
TARGET_LINK_LIBRARIES(libProjects libFORMAN)
TARGET_LINK_LIBRARIES(libProjects libUSAGE)
TARGET_LINK_LIBRARIES(libProjects libANN)
TARGET_LINK_LIBRARIES(libProjects libIOSUPPORT)

message(STATUS "Projects ends") #
