//
// Created by alex on 1/8/25.
//


#include "xx_tts.h"
#include "cgal_fixed_alphashape.h"
#include "TopoSegment.h"
#include <filesystem>

//#include <cstdio>

// get directory
string get_directory(const string &in_path) {
    // get absolute path
    std::filesystem::path path = std::filesystem::absolute(in_path);

    // return the parent path as a string
    return path.parent_path().string();
}

// get filename
string get_filename(const string &in_path) {
    // create a path object
    std::filesystem::path path(in_path);

    // return the filename as a string
    return path.filename().string();
}

int alpha_shape_generation(const string &infile, const double &alpha_sq_value, string &outfile) {
    string ab_filepath;
    std::filesystem::path absolutePath = std::filesystem::absolute(infile);
    ab_filepath = absolutePath.string();
    cout << "pts file path: " << ab_filepath << endl;
    cout << "alpha sq: " << alpha_sq_value << endl;
    // todo: update alpha shape code: better design and update the output file names
    string out_tocfile = "";
    std::size_t found2 = ab_filepath.find_last_of('.');
    // use the absolute path
    string outdir = ab_filepath.substr(0, found2); //
    stringstream ss;
    ss << fixed << setprecision(3) << outdir << "_a" << alpha_sq_value;
    if (ab_filepath.find(".pts") != string::npos) {
        ss << ".off";
    } else if (ab_filepath.find(".ply") != string::npos) {
        ss << ".ply";
    } else {
        cout << "no supported file format\n";
        return 0;
    }
    outfile = ss.str();
    cout << "output as3d: " << outfile << endl;
    //generate_fixed_alpha_shape(ab_filepath, alpha_sq_value, outfile, 3, out_tocfile);
    generate_fixed_alpha_shape_only(ab_filepath, alpha_sq_value, outfile);

    return 1;
}

int alpha_shape_segment(const string &infile,
                        const double &alpha_sq_value,
                        string &out_tocfile) {
    string ab_filepath;
    std::filesystem::path absolutePath = std::filesystem::absolute(infile);
    ab_filepath = absolutePath.string();
    cout << "pts file path: " << ab_filepath << endl;
    cout << "alpha sq: " << alpha_sq_value << endl;
    // if (outfile.empty()) {
    //     cout << "get components only\n";
    // } else {
    //     cout << "generate a complete alpha shape, too.\n";
    // }

    //string output_as_tocfile;
    // todo: update alpha shape code: better design and update the output file names
    generate_fixed_alpha_shape(ab_filepath, alpha_sq_value, "", 3, out_tocfile);
    std::cout << "out_tocfile: " << out_tocfile << endl;
    return 1;
}

int get_oversegments(const string &infile) {
    /*
    * input: alpa shape component alpha shape (i.e., .off) file
    * input parameter:
    */

    string workspace = get_directory(infile);
    string filename = get_filename(infile);
    cout << "workspace: " << workspace << endl;
    cout << "filename: " << filename << endl;

    IO_Timer timer;
    timer.start();
    // get over-segmentation results
    TopoSegment ts(infile, 3, false);

    string overseg_file = ""; // "_lbl.pts" // workspace + "/" + filename + "_oversegs.pts"; // over-segmentation xyzl file
    string forman_mins_file = ""; //workspace + "/" + filename + "_mins.off"; // forman mins file
    ts.cluster(overseg_file, forman_mins_file);

    timer.stop();
    cout << "\nover-segmentation total time: " << timer.getElapsedTimeInSec() << " s\n";

    return 0;
    
}

int extract_single_trees(const string &vg_mesh_file, const string &loc_file, string &outfile,
                         const double &th_p2trunk_distance,
                         const double &th_search_radius) {
    // infile: alpha shape file

    // extract single trees
    TopoSegment ts(vg_mesh_file, 3, false); // todo: use debug switch from the user
    // read trunk pts: xyzl.
    cout << "\nread seeds file\n";
    std::vector<std::vector<float> > detected_trunks;
    ts.tts_read_pts_fs(loc_file, detected_trunks);
    cout << "\ninit min-seeds\n";
    ts.initSeeds_mins_by_trunkpts_xyz(detected_trunks, th_p2trunk_distance, th_search_radius);
    cout << "\ngrow min-seeds\n";
    IO_Timer timer;
    timer.start();
    ts.growFromSeeds_mins_basic();
    timer.stop();
    cout << "grow time: " << timer.getElapsedTimeInSec() << " s\n";

    cout << "\nlabel pts\n";
    timer.start();
    //outfile = ts.label_pts_from_mins(true); // todo: add write file function.
    outfile = ts.label_pts_from_mins_parallel(true);
    timer.stop();
    cout << "label time: " << timer.getElapsedTimeInSec() << " s\n";

    return 1;
}
