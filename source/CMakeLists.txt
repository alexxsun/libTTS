cmake_minimum_required(VERSION 3.9)

add_subdirectory(iastar) # data structure encodes simplicial complex
add_subdirectory(forman) # core computation parts: Forman vpath, critical cells, etc
add_subdirectory(projects) # applications
add_subdirectory(usage)
add_subdirectory(ANN)

#add_subdirectory(cgal_as)