//
// Created by alex on 1/8/25.
//

#ifndef XX_TTS_H
#define XX_TTS_H

#include <string>
#include <vector>
using namespace std;

// generate alpha shape, no find component (no segmentation)
int alpha_shape_generation(const string &infile, const double &alpha_sq_value, string &outfile);

// segment entire points into subregions
// input: input pts file,
// input parameter: alpha value
// output: table of content file (file list file) of subregion pts file
int alpha_shape_segment(const string &infile, const double &alpha_sq_value, string &out_tocfile);

// wood-leaf separation
// input: alpa shape component alpha shape (i.e., .off) file
// input parameter:
// output: woody pts file
int wood_leaf_separation(const string &infile,
                         const double &th_knn_distance = 0.02,
                         string outfile = "");


// tree location detection
// input: input pts file
// input parameter:
// output: tree location file, data structure
// todo: implement here or in Python?
int detect_single_trees(const string &infile, const string &outfile,
                        std::vector<std::vector<float> > &detected_trunks);

// segment vegetation points to individual tree point clouds
// input: vegetation pts file, tree location file
// input parameter:
// output: labeled vg_file, xyzl
int extract_single_trees(const string &vg_mesh_file, const string &loc_file, string &outfile,
                         const double &th_p2trunk_distance = 0.2,
                         const double &th_search_radius = 0.25);

// int label_veg_points(const string &asfolder);


#endif //XX_TTS_H
