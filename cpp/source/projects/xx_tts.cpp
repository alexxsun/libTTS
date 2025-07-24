//
// Created by alex on 1/8/25.
//

#include <filesystem>
#include "xx_tts.h"
#include "cgal_fixed_alphashape_v2.h"
#include "TopoSegment.h"

int alpha_shape_generation(const string& infile, const double& alpha_sq_value, string& outfile) {
    string format = "off";

    if (outfile.empty()) {
        if (infile.find(".ply") != string::npos) {
            format = "ply";
        }
        std::stringstream ss;
        ss << std::filesystem::path(infile).stem().string()
            << "_a" << std::fixed << std::setprecision(3) << alpha_sq_value << "." << format;
        outfile = (std::filesystem::path(infile).parent_path() / ss.str()).string();
    } else {
        if (outfile.find(".ply") != string::npos) {
            format = "ply";
        }
    }
    generate_fixed_alpha_shape_only_v2(infile, alpha_sq_value, outfile, format);

    return 1;
}

int alpha_shape_segment(const string& infile,
                        const double& alpha_sq_value,
                        string& outdir) {
    string format = "off";
    if (infile.find(".ply") != string::npos) {
        format = "ply";
    }
    std::stringstream ss;
    ss << std::filesystem::path(infile).stem().string()
        << "_a" << std::fixed << std::setprecision(3) << alpha_sq_value << "." << format;
    //string outfile = (std::filesystem::path(infile).parent_path() / ss.str()).string();
    string outfile = ""; // we don't output entire alpha shape here

    ss.str("");
    ss.clear();
    ss << std::filesystem::path(infile).stem().string()
        << "_a" << std::fixed << std::setprecision(3) << alpha_sq_value << "_components";
    if (outdir.empty()) {
        outdir = (std::filesystem::path(infile).parent_path() / ss.str()).string();
    }
    std::cout << "out dir: " << outdir << "\n";
    // create outdir
    std::filesystem::path dir_path(outdir);
    if (std::filesystem::exists(dir_path)) {
        std::cout << "Directory '" << dir_path << "' already exists. Removing it." << std::endl;
        std::filesystem::remove_all(dir_path);
    }
    std::cout << "Creating directory '" << dir_path << "'." << std::endl;
    std::filesystem::create_directories(dir_path);

    //
    generate_fixed_alpha_shape_v2(infile, alpha_sq_value, outfile, outdir, format);
    return 1;
}

int get_oversegments(const string& infile, string& out_segfile) {
    /*
    * input: alpa shape component alpha shape (i.e., .off) file
    */

    IO_Timer timer;
    timer.start();
    // get over-segmentation results
    TopoSegment ts(infile, 3, false);

    // out_segfile: xx_lbl.pts/.ply
    ts.cluster(out_segfile);

    timer.stop();
    cout << "\nover-segmentation total time: " << timer.getElapsedTimeInSec() << " s\n";

    return 1;
}

int extract_single_trees(const string& vg_mesh_file, const string& loc_file, string& outfile,
                         const double& th_p2trunk_distance,
                         const double& th_search_radius) {
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
    //outfile = ts.label_pts_from_mins(true);
    outfile = ts.label_pts_from_mins_parallel(true);
    timer.stop();
    cout << "label time: " << timer.getElapsedTimeInSec() << " s\n";

    return 1;
}