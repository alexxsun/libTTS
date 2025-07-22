//
// Created by alex on 1/8/25.
//

#ifndef XX_TTS_H
#define XX_TTS_H

#include <string>
#include <vector>
using namespace std;

// generate alpha shape, no need to find components (no segmentation)
int alpha_shape_generation(const string& infile, const double& alpha_sq_value, string& outfile);

// segment entire points into subregions
int alpha_shape_segment(const string& infile, const double& alpha_sq_value, string& outdir);

// over-segmentation of alpha shape
int get_oversegments(const string& infile, string& out_segfile);

// segment vegetation points to individual tree point clouds
// input: vegetation pts file, tree location file
// input parameter:
// output: labeled vg_file, xyzl
int extract_single_trees(const string& vg_mesh_file, const string& loc_file, string& outfile,
                         const double& th_p2trunk_distance = 0.2,
                         const double& th_search_radius = 0.25);


#endif //XX_TTS_H