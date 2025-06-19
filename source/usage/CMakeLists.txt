cmake_minimum_required(VERSION 3.9)
message(STATUS "Usage starts...") #


find_package(Boost REQUIRED) # Boost not BOOST in Linux
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Boost_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Boost_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Boost_EXE_LINKER_FLAGS}")
endif ()



file(GLOB SOURCE_USAGE "*.h" "*.cpp")
add_library(libUSAGE STATIC ${SOURCE_USAGE})
message(STATUS "libUSAGE ends...") #

