cmake_minimum_required(VERSION 3.9)
message(STATUS "IASTAR starts...") #
#find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#    message(STATUS "OpenMP version: ${OpenMP_VERSION}")
#endif ()

find_package(Boost REQUIRED) # Boost not BOOST in Linux
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Boost_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Boost_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Boost_EXE_LINKER_FLAGS}")
endif ()

#find_package(OpenMP REQUIRED) # need it?
#if(OpenMP_CXX_FOUND)
#    message(STATUS "OpenMP version: ${OpenMP_VERSION}")
#endif()


file(GLOB SOURCE_IASTAR "*.h" "*.cpp")
add_library(libIASTAR STATIC ${SOURCE_IASTAR})
#target_link_libraries(libIASTAR PUBLIC OpenMP::OpenMP_CXX)
message(STATUS "IASTAR ends...") #

