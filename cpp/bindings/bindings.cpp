//
// Created by alex on 5/29/25.
//


#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For automatic conversion of std::string, std::vector, etc.
#include "../source/projects/xx_tts.h"      // TLS segmentation
#include "../source/projects/als_seg.h"  // ALS segmentation

namespace py = pybind11;

std::string generate_alpha_shape_py(
    const std::string& infile,
    double alpha_sq_value) {
    std::string output_filename_from_cpp;
    int status = alpha_shape_generation(infile, alpha_sq_value, output_filename_from_cpp);
    // Check the status code. If it's not success (1), raise an exception.
    if (status != 1) {
        throw std::runtime_error(
            "The C++ alpha_shape_generation function failed with status code " + std::to_string(status));
    }
    // On success, return the filename that the C++ function provided.
    return output_filename_from_cpp;
}

std::string get_oversegments_py(const string& infile) {
    std::string out_segfile;
    int status = get_oversegments(infile, out_segfile);
    if (status != 1) {
        throw std::runtime_error("The C++ function failed");
    }
    return out_segfile;
}

std::string tls_extract_single_trees_py(
    const std::string& vg_mesh_file,
    const std::string& loc_file,
    double th_p2trunk_distance,
    double th_search_radius) {
    std::string output_filename_from_cpp;
    int status = extract_single_trees(vg_mesh_file, loc_file, output_filename_from_cpp,
                                      th_p2trunk_distance, th_search_radius);

    if (status != 1) {
        throw std::runtime_error(
            "The C++ extract_single_trees function failed with status code " + std::to_string(status));
    }
    return output_filename_from_cpp;
}

PYBIND11_MODULE(_libtts, m) {
    // Python module name
    m.doc() = "Internal C++ bindings for libTTS"; // Optional module docstring

    // (infile, th_sq_value) -> outfile
    m.def("generate_alpha_shape_cpp",
          &generate_alpha_shape_py,
          R"doc(
            Generates a 3D alpha shape mesh from a point cloud.
            
            This function reads a .pts file, computes the alpha shape and saves the result to a new file.
            The alpha value controls the level of detail and concavity.
            
            Args:
                infile (str): Path to the input point cloud file (.pts).
                alpha_sq_value (float): The squared alpha value for the computation.
                    A smaller value results in a tighter mesh with more detail and
                    potential holes. A larger value results in a more convex shape.
                    Defaults to 0.01.
            
            Returns:
                str: The path to the generated output mesh file.
            )doc",
          py::arg("infile"),
          py::arg("alpha_sq_value") = 0.01
        );

    m.def("get_oversegments_cpp",
          &get_oversegments_py,
          R"doc(
            Performs over-segmentation on a 3D mesh file.
            
            This function takes a mesh file (e.g., from alpha shape generation)
            and divides it into a set of smaller subset made by points. 
            The result is a point cloud where each point is assigned a segment ID.
            
            Args:
                infile (str): Path to the input mesh file (.ply or .off).
            
            Returns:
                str: The path to the output file containing the segmented point
                    cloud with labels.
            )doc",
          py::arg("infile")
        );

    m.def("tls_extract_single_trees_cpp",
          &tls_extract_single_trees_py,
          R"doc(
            Extracts individual tree point clouds from the vegetation mesh.
            
            Using a vegetation mesh and a file of initial tree trunk locations, this
            function segments the mesh to isolate the points belonging to each
            individual tree.
            
            Args:
                vg_mesh_file (str): Path to the input vegetation mesh file.
                loc_file (str): Path to a file containing the X,Y,Z,L locations of
                    the tree trunks to be segmented.
                th_p2trunk_distance (float): The maximum distance a point can be
                    from a trunk to be considered part of that tree. Defaults to 0.2.
                th_search_radius (float): The search radius used during the
                    segmentation process. Defaults to 0.25.
            
            Returns:
                str: The path to the output file containing the points for all
                    extracted single trees, likely with labels.
            )doc",
          py::arg("vg_mesh_file"),
          py::arg("loc_file"),
          py::arg("th_p2trunk_distance") = 0.2,
          py::arg("th_search_radius") = 0.25
        );

    // ALS: segmentation. to improve.
    m.def("als_segment",
          &als_segment,
          R"doc(
            Segments Airborne Laser Scanning (ALS) data into individual tree clusters.
            
            This function processes an ALS point cloud, applying a segmentation
            algorithm based on persistence values to identify and separate
            individual tree crowns. It offers several thresholds to fine-tune the
            segmentation for different forest types.
            
            Args:
                infile (str): Path to the input ALS point cloud file.
                th_radius (float): Search radius parameter. Defaults to -1 (auto).
                th_forest_pers (float): Forest persistence threshold. Defaults to -1 (auto).
                th_height (float): Height threshold. Defaults to -1 (auto).
                th_pers_H (float): Persistence H threshold. Defaults to -1 (auto).
                th_pers_I (float): Persistence I threshold. Defaults to -1 (auto).
                th_height_funtype (bool): Flag for height function type. Defaults to True.
                ds (str): Optional output path for a downsampled cloud. Defaults to "".
            
            Returns:
                str: The path to the output file containing the segmented tree clusters.
            )doc",
          py::arg("infile"),
          py::arg("th_radius") = -1,
          py::arg("th_forest_pers") = -1,
          py::arg("th_height") = -1,
          py::arg("th_pers_H") = -1,
          py::arg("th_pers_I") = -1,
          py::arg("th_height_funtype") = true,
          py::arg("ds") = ""
        );
}