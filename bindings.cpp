//
// Created by alex on 5/29/25.
//


#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For automatic conversion of std::string, std::vector, etc.
#include "source/projects/xx_tts.h"      // TLS segmentation
#include "source/projects/als_seg.h"  // ALS segmentation

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

PYBIND11_MODULE(xx_tts_py, m) {
    // "xx_tts_py" will be the name of the Python module
    m.doc() = "Python bindings for the xx_tts C++ library"; // Optional module docstring

    // TLS: int alpha_shape_generation(const string &infile, const double &alpha_sq_value, string &outfile);
    m.def("alpha_shape_generation",
          &alpha_shape_generation_py,
          "Generates alpha shape. Returns a tuple (status_code, output_filename).",
          py::arg("infile"),
          py::arg("alpha_sq_value")
    );

    // TLS: int wood_leaf_separation(const string &infile, const double &th_knn_distance = 0.02, string outfile = "");
    m.def("wood_leaf_separation",
          &wood_leaf_separation,
          "Performs wood-leaf separation.",
          py::arg("infile"),
          py::arg("th_knn_distance") = 0.02, // Expose default value to Python
          py::arg("outfile") = "" // Expose default value to Python
    );

    // ALS segmentation
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
