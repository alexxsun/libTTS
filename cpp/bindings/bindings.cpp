//
// Created by alex on 5/29/25.
//


#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For automatic conversion of std::string, std::vector, etc.
#include "../source/projects/xx_tts.h"      // TLS segmentation
#include "../source/projects/als_seg.h"  // ALS segmentation

namespace py = pybind11;

// Wrapper for alpha_shape_generation
// The C++ function `alpha_shape_generation` takes `outfile` as a reference
// parameter to return a value. In Python, it's more idiomatic for functions
// to return values directly. This wrapper facilitates that.
// It will return a pair: (int_status_code, output_filename_string)
std::pair<int, std::string> alpha_shape_generation_py(const std::string &infile, const double &alpha_sq_value) {
    std::string outfile_str; // Variable to capture the output filename
    int result = alpha_shape_generation(infile, alpha_sq_value, outfile_str);
    return std::make_pair(result, outfile_str);
}

PYBIND11_MODULE(libtts_py, m) {
    // Python module name: libtts_py
    m.doc() = "Python bindings for the xx_tts C++ library"; // Optional module docstring

    // TLS: int alpha_shape_generation(const string &infile, const double &alpha_sq_value, string &outfile);
    m.def("alpha_shape_generation",
          &alpha_shape_generation_py,
          "Generates 3D alpha shape. Returns a tuple (status_code, output_filename).",
          py::arg("infile"),
          py::arg("alpha_sq_value")
    );

    // TLS: int get_oversegments(const string &infile);
    // todo: update the function ==> TLS: over-segmentation
    m.def("get_oversegments",
          &get_oversegments,
          "Performs over-segmentation.",
          py::arg("infile")
    );

    // TLS: int extract_single_trees(const string &vg_mesh_file, const string &loc_file, string &outfile,
    //                     const double &th_p2trunk_distance = 0.2,
    //                     const double &th_search_radius = 0.25);
    m.def("extract_single_trees",
          &extract_single_trees,
          "Segments vegetation points into individual tree point clouds.",
          py::arg("vg_mesh_file"),
          py::arg("loc_file"),
          py::arg("outfile"),
          py::arg("th_p2trunk_distance") = 0.2,
          py::arg("th_search_radius") = 0.25
    );

    // ALS: segmentation
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
