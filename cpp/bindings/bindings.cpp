//
// Created by alex on 5/29/25.
//


#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For automatic conversion of std::string, std::vector, etc.
#include "../source/projects/xx_tts.h"      // TLS segmentation
#include "../source/projects/als_seg.h"  // ALS segmentation

namespace py = pybind11;

std::string generate_alpha_shape_py(
    const std::string &infile,
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

std::string get_oversegments_py(const string &infile) {
    std::string out_segfile, out_minsfile;
    int status = get_oversegments(infile, out_segfile, out_minsfile);
    if (status != 1) {
        throw std::runtime_error("The C++ function failed");
    }
    return out_segfile;
}

std::string tls_extract_single_trees_py(
    const std::string &vg_mesh_file,
    const std::string &loc_file,
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
          "Generates 3D alpha shape. infile: pts file",
          py::arg("infile"),
          py::arg("alpha_sq_value") = 0.01
    );


    m.def("get_oversegments_cpp",
          &get_oversegments_py,
          "Performs over-segmentation. infile: shape file",
          py::arg("infile")
    );


    m.def("tls_extract_single_trees_cpp",
          &tls_extract_single_trees_py,
          "Segments trees from TLS data. Infile: vg_mesh_file",
          py::arg("vg_mesh_file"),
          py::arg("loc_file"),
          py::arg("th_p2trunk_distance") = 0.2,
          py::arg("th_search_radius") = 0.25
    );

    // ALS: segmentation. to improve.
    m.def("als_segment",
          &als_segment,
          "Segments ALS data into tree clusters based on persistence values.",
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
