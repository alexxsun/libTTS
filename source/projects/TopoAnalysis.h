//
// Created by alex on 1/2/24.
//

#ifndef TTS_TOPOANALYSIS_H
#define TTS_TOPOANALYSIS_H

#include "../forman/formangradient.h"
#include <filesystem> // C++ 17

// : public FormanGradient
class TopoAnalysis {

protected:
    bool _showlog = false;
    string _infile;
    string _file_name;
    string _workspace;

    void _set_workspace(const string &infile) {
        _infile = infile;
        std::filesystem::path absolutePath = std::filesystem::absolute(infile);
        _file_name = absolutePath.stem().string();
        std::filesystem::path parentPath = absolutePath.parent_path();
        _workspace = parentPath.string() + "/";

        cout << "infile: " << _infile << endl;
        cout << "file name: " << _file_name << endl;
        cout << "workspace: " << _workspace << endl;

    }


public:
    // : FormanGradient(infile, fid)
    TopoAnalysis(const string &infile, const int &fid, const bool &in_debug = false) {

        // not use the FormanGradient as the base class, because
        // it will first finish everything in FormanGradient.
        // May need to update the FormanGradient class later.
        _showlog = in_debug;
        _set_workspace(infile);

//        if (sc.getVertexCoordSize() != 3) {
//            std::cout << "please ensure input point: x y z\n";
//            exit(1); // sth wrong
//        }

        // let's have some basic info based on Forman Gradient
        // to check Ciccio's repo code.
        // https://github.com/IuricichF/SimplForman3D

        FormanGradient fg(infile, fid);

        // compute the Forman gradient
        std::cout << "compute forman gradient\n";
        fg.computeFormanGradient(true);

        std::cout << "Output Forman Gradient related results\n";
        // input mesh
        fg.saveScalarFieldVTK_OG((_workspace + _file_name + "_scalarField.vtk").c_str());
        // critical cells
        fg.xx_output_critical_cells(_workspace + _file_name + "_critical_cells.pts");
        fg.xx_output_critical_cells_vtk(_workspace + _file_name + "_critical_cells.vtk");

        // critical cell persistence. Note: not correct.
        // fg.write_critical_persistence_2d(_workspace + _file_name + "_critical_pers.txt");

        // v-path
        // fg.visMorse(); // descending1manifold, ascending1cells_correct.vtk
        fg.xx_visMorse(_workspace + _file_name); // descending1manifold, ascending1cells_correct.vtk
        fg.xx_vis_VPath_01(_workspace + _file_name + "_vpaths_01.vtk", false, false);
        fg.xx_vis_VPath_12(_workspace + _file_name + "_vpaths_12.vtk", false, false);

        // network of critical cells
        fg.xx_critical_cells_net_2d(_workspace + _file_name + "_cp_net.vtk", _workspace + _file_name + "_cp_net.txt");

        // gradient. note: nor correct. not ready.
        //fg.writeVTK_gradient_2d(_workspace + _file_name + "_gradient_2d.vtk");

        // old code...

        fg.xx_vis_CriticalCells_01(_workspace + _file_name + "_critical_cells_01.vtk", false);


    }

};


#endif //TTS_TOPOANALYSIS_H
