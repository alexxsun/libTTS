//
// Created by alex on 12/26/23.
//

#include "TopoSegment.h"
#include "../io_support/happly.h"


void TopoSegment::_get_connected_points(const double &th_h1, const double &th_h2, vector<std::set<int> > &con_pts) {
    /* generate vv relationship of points (x,y,z,...) with p.z in [th_h1, th_h2]
        * point: xyzv
        *
        */
    if (_showlog) cout << "get connected low points in [" << th_h1 << ", " << th_h2 << "]\n";

    std::map<int, set<int> > vv; //
    std::map<int, int> to_process; // 1: processed, otherwise 0
    for (int i = 0; i < sc.getVerticesNum(); ++i) {
        Vertex v = sc.getVertex(i);
        if (v.getCoordinate(2) < th_h1 || v.getCoordinate(2) > th_h2) {
            // xyz
            continue;
        }
        to_process[i] = 0;
        implicitS tmp(i);
        for (auto nbh: *sc.adjacents(tmp)) {
            int j = nbh.getConstVertices()[0];
            Vertex v2 = sc.getVertex(j);
            if (v2.getCoordinate(2) < th_h1 || v2.getCoordinate(2) > th_h2) {
                // xyz
                continue;
            }
            if (i == j) {
                continue;
            }
            vv[i].insert(j);
            to_process[j] = 0;
        }
    }
    if (_showlog) cout << "vv done\n";

    // DFS?
    for (const auto &pid: to_process) {
        if (pid.second > 0) {
            // processed
            continue;
        }
        // a new set starts
        int vid = pid.first;
        to_process[vid] = 1;
        std::set<int> tmp_con = {vid}; // init connected points
        std::vector<int> nbhs; // init neighbors.
        for (const auto &nvid: vv[vid]) {
            nbhs.push_back(nvid);
        }
        while (!nbhs.empty()) {
            int nbh = nbhs.back();
            nbhs.pop_back();
            if (to_process[nbh] > 0) {
                // processed
                continue;
            }
            to_process[nbh] = 1;
            tmp_con.insert(nbh);
            for (const auto &nvid: vv[nbh]) {
                nbhs.push_back(nvid);
            }
        }
        //
        con_pts.push_back(tmp_con);
    }

    if (_showlog) cout << "connected low points done\n";
}


std::vector<int> TopoSegment::_label_points_by_grouped_mins(const map<int, std::vector<implicitS> > &gp_mins) {
    /*
    * label: mins
    * -2: boundary, 0: unlabeled, >=1: valid labels
    */
    if (_showlog) cout << "\nlabel points\n";

    std::vector<int> pts_lbls(sc.getVerticesNum(), 0);

    for (const auto &gps: gp_mins) {
        int lbl = gps.first;
        if (lbl < 1) {
            continue; // valid lbl starts from 1
        }
        for (const auto &cc: gps.second) {
            SSet tmp_cells;
            computeAscendingCell(true, cc, tmp_cells);
            // label points in the ascending/descending cells
            for (const auto &ac: tmp_cells) {
                for (const auto &vid: ac.getConstVertices()) {
                    if (pts_lbls[vid] > 0) {
                        int old_lbl = pts_lbls[vid];
                        if (old_lbl != lbl) {
                            pts_lbls[vid] = -2; //  point is labeled at least twice
                        }
                    } else {
                        pts_lbls[vid] = lbl;
                    }
                }
            }
        }
    } // label points done.

    if (_showlog) {
        cout << "# of final labels: " << gp_mins.size() << endl;
        std::vector<std::pair<int, int> > lblnums;
        for (const auto &gps: gp_mins) {
            int n = count(pts_lbls.begin(), pts_lbls.end(), gps.first);
            // cout << "lbl: " << gps.first << ", # " << n << endl;
            lblnums.push_back(std::pair<int, int>(gps.first, n));
        }
        std::sort(lblnums.begin(), lblnums.end(),
                  [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.second > b.second; });
        for (const auto &ln: lblnums) {
            cout << "lbl: " << ln.first << ", # " << ln.second << endl;
        }
    }

    return pts_lbls;
}

int TopoSegment::_label_points_by_grouped_mins_parallel(const map<int, std::vector<implicitS> > &gp_mins) {
    /*
    * label: mins
    * -2: boundary, 0: unlabeled, >=1: valid labels
    */
    if (_showlog) cout << "\nlabel points in parallel\n";

    _pts_lbls = std::vector<int>(sc.getVerticesNum(), 0);

    // Convert map to vector for better parallel processing
    std::vector<std::pair<int, std::vector<implicitS> > > gp_mins_vec(gp_mins.begin(), gp_mins.end());

#pragma omp parallel num_threads(6)
    {
        // Create thread-local vector for labels
        std::vector<int> local_labels(sc.getVerticesNum(), 0);

#pragma omp for schedule(dynamic) nowait
        for (size_t i = 0; i < gp_mins_vec.size(); ++i) {
            int lbl = gp_mins_vec[i].first;
            if (lbl < 1) continue; // valid lbl starts from 1

            for (const auto &cc: gp_mins_vec[i].second) {
                SSet tmp_cells;
                computeAscendingCell(true, cc, tmp_cells);

                // Label points in the ascending/descending cells
                for (const auto &ac: tmp_cells) {
                    for (const auto &vid: ac.getConstVertices()) {
                        if (local_labels[vid] > 0) {
                            int old_lbl = local_labels[vid];
                            if (old_lbl != lbl) {
                                local_labels[vid] = -2; // point is labeled at least twice
                            }
                        } else {
                            local_labels[vid] = lbl;
                        }
                    }
                }
            }
        }

        // Merge thread-local results
#pragma omp critical
        {
            for (size_t i = 0; i < _pts_lbls.size(); ++i) {
                if (local_labels[i] != 0) {
                    if (_pts_lbls[i] > 0) {
                        if (_pts_lbls[i] != local_labels[i]) {
                            _pts_lbls[i] = -2;
                        }
                    } else {
                        _pts_lbls[i] = local_labels[i];
                    }
                }
            }
        }
    }


    if (_showlog && false) {
        // skip for the test
        cout << "# of final labels: " << gp_mins.size() << endl;
        std::vector<std::pair<int, int> > lblnums;
        for (const auto &gps: gp_mins) {
            int n = count(_pts_lbls.begin(), _pts_lbls.end(), gps.first);
            // cout << "lbl: " << gps.first << ", # " << n << endl;
            lblnums.push_back(std::pair<int, int>(gps.first, n));
        }
        std::sort(lblnums.begin(), lblnums.end(),
                  [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.second > b.second; });
        for (const auto &ln: lblnums) {
            cout << "lbl: " << ln.first << ", # " << ln.second << endl;
        }
    }

    return 1;
}

void TopoSegment::_output_pts_with_label_pts(const string &outfile, const vector<int> &lbls, const bool &scaled) {
    /*
        * x y z lbl
        * ...
        */
    if (_showlog) cout << "write: " << outfile << endl;

    ofstream ofs(outfile);
    ofs << fixed << setprecision(3);
    //  ofs << "# OFF\n";
    //  ofs << sc.getVerticesNum() << " 0 0\n";

    for (int i = 0; i < sc.getVerticesNum(); ++i) {
        Vertex v = sc.getVertex(i);
        if (scaled) {
            ofs << v.getCoordinate(0) - sc.sc_min_x << " " << v.getCoordinate(1) - sc.sc_min_y << " "
                    << v.getCoordinate(2)
                    << " " << lbls[i] << endl;
        } else {
            ofs << v.getCoordinate(0) << " " << v.getCoordinate(1) << " " << v.getCoordinate(2) << " " << lbls[i]
                    << endl;
        }
    }
}


void TopoSegment::_output_pts_with_label_pts_ply(const string &outfile, const vector<int> &lbls,
                                                 const bool &scaled) {
    /*
        * x y z lbl
        * ...
        */
    if (_showlog) cout << "write: " << outfile << endl;
    std::vector<std::array<double, 3> > ply_vertices;
    for (int i = 0; i < sc.getVerticesNum(); ++i) {
        Vertex v = sc.getVertex(i);
        if (scaled) {
            ply_vertices.push_back({
                v.getCoordinate(0) - sc.sc_min_x, v.getCoordinate(1) - sc.sc_min_y,
                v.getCoordinate(2) - sc.sc_min_z
            });
        } else {
            ply_vertices.push_back({v.getCoordinate(0), v.getCoordinate(1), v.getCoordinate(2)});
        }
    }

    happly::PLYData plyOut;
    plyOut.addVertexPositions(ply_vertices);
    // todo: make a good name, check the concept of element, property...
    plyOut.getElement("vertex").addProperty<int>("label", lbls);
    // Write binary
    plyOut.write(outfile, happly::DataFormat::Binary);
}

void TopoSegment::_output_pts_with_label_vtk(const string &outfile, const vector<int> &lbls, const bool &scaled) {
    ofstream ofs(outfile);
    ofs << fixed << setprecision(3);
    ofs << "# vtk DataFile Version 2.0\n"
            "\n"
            "ASCII \n"
            "DATASET UNSTRUCTURED_GRID\n"
            "\n";
    ofs << "POINTS " << sc.getVerticesNum() << " float\n";

    for (int i = 0; i < sc.getVerticesNum(); ++i) {
        Vertex v = sc.getVertex(i);
        if (scaled) {
            ofs << v.getCoordinate(0) - sc.sc_min_x << " " << v.getCoordinate(1) - sc.sc_min_y << " "
                    << v.getCoordinate(2)
                    << endl;
        } else {
            ofs << v.getCoordinate(0) << " " << v.getCoordinate(1) << " " << v.getCoordinate(2) << endl;
        }
    }

    ofs << "CELLS " << sc.getVerticesNum() << " " << 2 * sc.getVerticesNum() << endl;
    for (int i = 0; i < sc.getVerticesNum(); ++i) {
        ofs << "1 " << i << endl;
    }
    ofs << "CELL_TYPES " << sc.getVerticesNum() << endl;
    // VTK_VERTEX (=1)
    for (int i = 0; i < sc.getVerticesNum(); ++i) {
        ofs << "1 " << endl;
    }

    ofs << "CELL_DATA " << sc.getVerticesNum() << "\n";
    ofs << "FIELD FieldData 1\n";
    ofs << "lbl 1 " << sc.getVerticesNum() << " int" << endl;
    for (int i = 0; i < sc.getVerticesNum(); ++i) {
        ofs << lbls[i] << endl;
    }
}

void TopoSegment::_output_mins_pts(const std::vector<implicitS> &mins, const std::string &ptsfile, const bool &scaled) {
    ofstream ofs(ptsfile);
    ofs << fixed << setprecision(2);

    // points
    for (const auto &v: mins) {
        Vertex ver = sc.getVertex(v.getConstVertices()[0]);
        if (scaled) {
            ofs << ver.getCoordinate(0) - sc.sc_min_x << " " << ver.getCoordinate(1) - sc.sc_min_y << " "
                    << ver.getCoordinate(2) << endl;
        } else {
            ofs << ver.getCoordinate(0) << " " << ver.getCoordinate(1) << " " << ver.getCoordinate(2) << endl;
        }
    }
}

void TopoSegment::_output_mins_vtk(const std::vector<implicitS> &mins, const std::string &vtktfile,
                                   const bool &scaled) {
    int pts_num = mins.size();
    int cell_num = mins.size();
    ofstream ofs(vtktfile);
    ofs << fixed << setprecision(2);
    ofs << "# vtk DataFile Version 2.0\n"
            "tm graph\n"
            "ASCII\n"
            "DATASET UNSTRUCTURED_GRID\n";
    // points
    ofs << "POINTS " << pts_num << " float\n";
    for (const auto &v: mins) {
        Vertex ver = sc.getVertex(v.getConstVertices()[0]);
        if (scaled) {
            ofs << ver.getCoordinate(0) - sc.sc_min_x << " " << ver.getCoordinate(1) - sc.sc_min_y << " "
                    << ver.getCoordinate(2) << endl;
        } else {
            ofs << ver.getCoordinate(0) << " " << ver.getCoordinate(1) << " " << ver.getCoordinate(2) << endl;
        }
    }
    // cells: VTK_VERTEX (=1)
    ofs << "CELLS " << cell_num << " " << cell_num * 2 << endl;
    for (int i = 0; i < pts_num; ++i) {
        ofs << 1 << " " << i << endl;
    }

    ofs << "CELL_TYPES " << cell_num << endl;
    for (int i = 0; i < cell_num; ++i) {
        ofs << 1 << "\n";
    }
    ofs << endl;

    // field data
    ofs << "POINT_DATA " << cell_num << endl;
    ofs << "FIELD FieldData 2\n";
    ofs << "pid 1 " << cell_num << " int\n";
    for (int i = 0; i < cell_num; ++i) {
        ofs << i << " \n";
    }
    ofs << endl;
    ofs << "v 1 " << cell_num << " float\n";
    for (const auto &v: mins) {
        Vertex ver = sc.getVertex(v.getConstVertices()[0]);
        ofs << ver.getCoordinate(3) << endl;
    }

    ofs.close();
}

void TopoSegment::_output_gpmins_pts(const map<int, std::vector<implicitS> > &gp_mins, const std::string &ptsfile,
                                     const bool &scaled) {
    // gp_min: lbl: vector of mins
    ofstream ofs(ptsfile);
    ofs << fixed << setprecision(2);

    // points
    for (const auto &gps: gp_mins) {
        int lbl = gps.first;
        for (const auto &v: gps.second) {
            Vertex ver = sc.getVertex(v.getConstVertices()[0]);
            if (scaled) {
                ofs << ver.getCoordinate(0) - sc.sc_min_x << " " << ver.getCoordinate(1) - sc.sc_min_y << " "
                        << ver.getCoordinate(2) << " " << lbl << endl;
            } else {
                ofs << ver.getCoordinate(0) << " " << ver.getCoordinate(1) << " " << ver.getCoordinate(2) << " " << lbl
                        << endl;
            }
        }
    }
    ofs.close();
}
