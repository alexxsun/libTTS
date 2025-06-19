cmake_minimum_required(VERSION 3.1)
message(STATUS "ANN starts...") #

#file(GLOB SOURCE_IASTAR "*.h" "*.cpp")
include_directories(.)
#add_subdirectory(src)
#add_subdirectory(ANN)
set(SOURCE_ANN src/ANN.cpp
        src/bd_fix_rad_search.cpp
        src/bd_pr_search.cpp
        src/bd_search.cpp
        src/bd_tree.cpp
        src/bd_tree.h
        src/brute.cpp
        src/kd_dump.cpp
        src/kd_fix_rad_search.cpp
        src/kd_fix_rad_search.h
        src/kd_pr_search.cpp
        src/kd_pr_search.h
        src/kd_search.cpp
        src/kd_search.h
        src/kd_split.cpp
        src/kd_split.h
        src/kd_tree.cpp
        src/kd_tree.h
        src/kd_util.cpp
        src/kd_util.h
        src/perf.cpp
        src/pr_queue_k.h
        src/pr_queue.h
        ANN/ANN.h
        ANN/ANNperf.h
        ANN/ANNx.h)
add_library(libANN STATIC ${SOURCE_ANN})
message(STATUS "ANN ends...") #



