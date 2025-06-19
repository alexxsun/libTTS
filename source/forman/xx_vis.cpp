//
// Created by alex on 3/3/21.
// Visualization code by Xin.

//
#include "formangradient.h"

void FormanGradient::xx_output_critical_cells(const string &outfile) {
    // min: vertex, saddle: edge, middle point, triangle: centroid point
    // output format
    // x y z type [min:0, saddle:1, ...]
    std::vector<std::array<double, 4>> cc_pts;  // (x,y,z,type)
    for (int i = 0; i < 3; i++) {
        for (const auto &cc: criticalS[i]) {
            vector<int> vertices = cc.getConstVertices();
            //
            if (i == 0) {
                Vertex v = sc.getVertex(vertices[0]);
                //
                double x, y, z;
                x = v.getCoordinate(0);
                y = v.getCoordinate(1);
                z = v.getCoordinate(2);
                cc_pts.push_back({x, y, z, static_cast<double>(i)});
            }
            //
            if (i == 1) {
                Vertex v = sc.getVertex(vertices[0]);
                Vertex v2 = sc.getVertex(vertices[1]);
                v.middlePoint(v2);
                //
                double x, y, z;
                x = v.getCoordinate(0);
                y = v.getCoordinate(1);
                z = v.getCoordinate(2);
                cc_pts.push_back({x, y, z, static_cast<double>(i)});
            }
            if (i == 2) {
                Vertex v = sc.getVertex(vertices[0]);
                Vertex v2 = sc.getVertex(vertices[1]);
                Vertex v3 = sc.getVertex(vertices[2]);
                double x, y, z;
                x = (v.getCoordinate(0) + v2.getCoordinate(0) + v3.getCoordinate(0)) / 3.0;
                y = (v.getCoordinate(1) + v2.getCoordinate(1) + v3.getCoordinate(1)) / 3.0;
                z = (v.getCoordinate(2) + v2.getCoordinate(2) + v3.getCoordinate(2)) / 3.0;
                cc_pts.push_back({x, y, z, static_cast<double>(i)});
            }
        }
    }

    // write file
    ofstream ofs(outfile);
    ofs << fixed << setprecision(2);
    for (const auto &p: cc_pts) {
        ofs << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
    }
    ofs.close();
}

int FormanGradient::xx_output_critical_cells_vtk(const string &vtkfile) {
    // min: vertex, saddle: edge, max: triangle

    std::vector<std::array<double, 3>> cc_pts;
    int cp_min_num = 0, cp_min_st = 0, cp_saddle_num = 0, cp_saddle_st = 0, cp_max_num = 0, cp_max_st = 0;

    // mins
    cp_min_num = criticalS[0].size();
    for (const auto &cc: criticalS[0]) {
        vector<int> vertices = cc.getConstVertices();
        Vertex v = sc.getVertex(vertices[0]);
        double x, y, z;
        x = v.getCoordinate(0);
        y = v.getCoordinate(1);
        z = v.getCoordinate(2);
        cc_pts.push_back({x, y, z});
    }
    // saddles
    cp_saddle_st = cc_pts.size();
    cp_saddle_num = criticalS[1].size();
    for (const auto &cc: criticalS[1]) {
        vector<int> vertices = cc.getConstVertices();
        Vertex v = sc.getVertex(vertices[0]);
        Vertex v2 = sc.getVertex(vertices[1]);
        double x, y, z;
        x = v.getCoordinate(0);
        y = v.getCoordinate(1);
        z = v.getCoordinate(2);
        cc_pts.push_back({x, y, z});
        x = v2.getCoordinate(0);
        y = v2.getCoordinate(1);
        z = v2.getCoordinate(2);
        cc_pts.push_back({x, y, z});
    }
    // maxs
    cp_max_st = cc_pts.size();
    cp_max_num = criticalS[2].size();
    for (const auto &cc: criticalS[2]) {
        vector<int> vertices = cc.getConstVertices();
        Vertex v = sc.getVertex(vertices[0]);
        Vertex v2 = sc.getVertex(vertices[1]);
        Vertex v3 = sc.getVertex(vertices[2]);
        double x, y, z;
        x = v.getCoordinate(0);
        y = v.getCoordinate(1);
        z = v.getCoordinate(2);
        cc_pts.push_back({x, y, z});
        //
        x = v2.getCoordinate(0);
        y = v2.getCoordinate(1);
        z = v2.getCoordinate(2);
        cc_pts.push_back({x, y, z});
        //
        x = v3.getCoordinate(0);
        y = v3.getCoordinate(1);
        z = v3.getCoordinate(2);
        cc_pts.push_back({x, y, z});
    }

    int cells_num = cp_min_num + cp_saddle_num + cp_max_num;

    // write file
    ofstream ofs(vtkfile);
    ofs << fixed << setprecision(3);
    ofs << "# vtk DataFile Version 2.0\n"
           "critical cells\n"
           "ASCII\n"
           "DATASET UNSTRUCTURED_GRID\n";
    // points
    ofs << "POINTS " << cc_pts.size() << " float\n";
    for (const auto &pt: cc_pts) {
        ofs << pt[0] << " " << pt[1] << " " << pt[2] << endl;
    }
    // cells
    ofs << "CELLS " << cells_num << " "
        << 2 * cp_min_num + 3 * cp_saddle_num + 4 * cp_max_num
        << endl;
    for (int i = 0; i < cp_min_num; ++i) {
        ofs << "1 " << i << endl;
    }
    for (int i = 0; i < cp_saddle_num; ++i) {
        ofs << "2 " << cp_saddle_st + 2 * i << " " << cp_saddle_st + 2 * i + 1 << endl;
    }
    for (int i = 0; i < cp_max_num; ++i) {
        ofs << "3 " << cp_max_st + 3 * i << " " << cp_max_st + 3 * i + 1 << " " << cp_max_st + 3 * i + 2 << endl;
    }

    // VTK_VERTEX (=1), VTK_LINE (=3), VTK_TRIANGLE(=5)
    ofs << "CELL_TYPES " << cells_num << endl;
    for (auto i = 0; i < cp_min_num; ++i) {
        ofs << "1 ";
    }
    ofs << endl;
    for (auto i = 0; i < cp_saddle_num; ++i) {
        ofs << "3 ";
    }
    ofs << endl;
    for (int i = 0; i < cp_max_num; ++i) {
        ofs << "5 ";
    }
    ofs << endl;
    // field data: cell
    ofs << "CELL_DATA " << cells_num << endl;
    ofs << "FIELD FieldData 1\n";
    ofs << "cp_type 1 " << cells_num << " int\n";
    for (auto i = 0; i < cp_min_num; ++i) {
        ofs << "0 ";
    }
    ofs << endl;
    for (auto i = 0; i < cp_saddle_num; ++i) {
        ofs << "1 ";
    }
    ofs << endl;
    for (auto i = 0; i < cp_max_num; ++i) {
        ofs << "2 ";
    }
    ofs << endl;
    //
    ofs.close();

    return 0;// everything goes well
}

void FormanGradient::xx_vis_CriticalCells_01(const string &vtkfile, const bool &scaled) {
    // scaled: scaled version or not. only scale x,y
    //

    vector<vector<int>> ccs;       // critical cells
    vector<Vertex> saddle_points;  // saddle middle point

    for (const auto &cc: criticalS[0]) {
        ccs.push_back(cc.getConstVertices());
    }
    size_t saddle_sid = ccs.size();  // saddle begin id, also the pts number of the
    // critical minimal points
    for (const auto &cc: criticalS[1]) {
        ccs.push_back(cc.getConstVertices());
        vector<int> vertices = cc.getConstVertices();
        Vertex v = sc.getVertex(vertices[0]);
        Vertex v2 = sc.getVertex(vertices[1]);
        v.middlePoint(v2);
        saddle_points.push_back(v);
    }

    map<int, int> old2new;  // original point id -> written vtk point id
    int pts_num = 0;
    for (auto &cp: ccs) {
        for (auto p: cp) {
            if (old2new.find(p) == old2new.end()) {
                old2new[p] = pts_num;
                pts_num += 1;
            }
        }
    }
    vector<int> new2old(pts_num);
    for (auto &on: old2new) {
        new2old[on.second] = on.first;
    }

    size_t cell_num = saddle_sid + 2 * saddle_points.size();

    // count # of points below each critical point
    vector<int> lp_nums;  // low points number todo: remove it.
    for (int i = 0; i < saddle_sid; ++i) {
        int nid = ccs[i][0];
        Vertex p = sc.getVertex(nid);
        lp_nums.push_back(0);  // lp_nums.push_back(tm_getLPnumber(p));
    }

    for (auto &v: saddle_points) {
        lp_nums.push_back(0);  // lp_nums.push_back(tm_getLPnumber(v));
        lp_nums.push_back(0);  // lp_nums.push_back(tm_getLPnumber(v));
    }

    // write vtk
    ofstream ofs(vtkfile);
    ofs << "# vtk DataFile Version 2.0\n"
           "critical cells\n"
           "ASCII\n"
           "DATASET UNSTRUCTURED_GRID\n";
    // points
    ofs << "POINTS " << pts_num + saddle_points.size() << " float\n";
    for (const auto &v: new2old) {
        Vertex ver = sc.getVertex(v);
        if (scaled) {
            ofs << ver.getCoordinate(0) - sc.sc_min_x << " " << ver.getCoordinate(1) - sc.sc_min_y << " "
                << ver.getCoordinate(2) << endl;
        } else {
            ofs << ver.getCoordinate(0) << " " << ver.getCoordinate(1) << " " << ver.getCoordinate(2) << endl;
        }
    }
    for (auto &ver: saddle_points) {
        if (scaled) {
            ofs << ver.getCoordinate(0) - sc.sc_min_x << " " << ver.getCoordinate(1) - sc.sc_min_y << " "
                << ver.getCoordinate(2) << endl;
        } else {
            ofs << ver.getCoordinate(0) << " " << ver.getCoordinate(1) << " " << ver.getCoordinate(2) << endl;
        }
    }

    // cells: (two points) VTK_LINE = 3, VTK_VERTEX =1
    ofs << "CELLS " << cell_num << " " << 2 * (saddle_sid + saddle_points.size()) + 3 * saddle_points.size() << endl;
    for (auto i = 0; i < saddle_sid; ++i) {
        ofs << "1 ";
        for (auto &v: ccs[i]) {
            ofs << old2new[v] << " ";
        }
        ofs << endl;
    }
    int saddle_sid2 = pts_num;
    for (auto i = 0; i < saddle_points.size(); ++i) {
        ofs << "1 " << saddle_sid2 << endl;
        saddle_sid2 += 1;
    }
    for (size_t i = saddle_sid; i < ccs.size(); ++i) {
        ofs << "2 ";
        for (auto &v: ccs[i]) {
            ofs << old2new[v] << " ";
        }
        ofs << endl;
    }

    ofs << "CELL_TYPES " << cell_num << endl;
    for (auto i = 0; i < saddle_sid; ++i) {
        ofs << "1 ";
    }
    ofs << endl;
    for (auto i = 0; i < saddle_points.size(); ++i) {
        ofs << "1 ";
    }
    ofs << endl;
    for (int i = 0; i < saddle_points.size(); ++i) {
        ofs << "3 ";
    }
    ofs << endl;

    // field data: point - need it?
    //    ofs << "POINT_DATA " << saddle_sid + saddle_points.size() << endl;;
    //    ofs << "FIELD FieldData 1\n";
    //    ofs << "cp_type 1 " << saddle_sid + saddle_points.size() << " float\n";
    //    for (size_t i = 0; i < saddle_sid; ++i) {
    //        ofs << "0 ";
    //    }
    //    for (int i = 0; i < saddle_points.size(); ++i) {
    //        ofs << "1 ";
    //    }
    //    ofs << endl;

    // field data: cell
    ofs << "CELL_DATA " << cell_num << endl;
    ofs << "FIELD FieldData 2\n";
    ofs << "cp_type 1 " << cell_num << " int\n";
    for (auto i = 0; i < saddle_sid; ++i) {
        ofs << "0 ";
    }
    for (auto i = 0; i < 2 * saddle_points.size(); ++i) {
        ofs << "1 ";
    }
    ofs << endl;
    ofs << "num 1 " << cell_num << " int\n";
    for (auto &n: lp_nums) {
        ofs << n << " ";
    }
    ofs << endl;
    ofs.close();
}

void FormanGradient::xx_vis_paired_arrows(const string &outfile) {
    // todo: to finish it

    // generate arrow csv file
    // x,y,z,dx,dy,dz (x,y,y): start point, (dx,dy,dz): arrow direction
    // how to use?
    // read arrow csv in paraview -> table2point -> calculator (calculate direction vector from (dx,dy,dz)) -> Glyph

    // get info

    // write
}

void FormanGradient::xx_vis_VPath_01(const string &vtkfile, const bool &scaled, const bool &outdebugfile) {
    // get separatrix v-path: a array of cells. -> save vertices of cells
    // -> vector<vector<int>> vpaths; start from the end point of saddle edge
    vector<vector<int>> v_paths;
    for (const auto &cc_1: criticalS[1]) {  // critical 1-cells
        for (int i = 0; i < 2; ++i) {          // cell.getConstVertices().size()
            vector<int> vpath;
            int cv = cc_1.getConstVertices()[i];
            vpath.push_back(cv);

            implicitS s(cv);  // vertex simplex
            implicitS next;
            // start from its vertex
            if (getPair(s, next)) {
                // todo: do we need to "assert" here? It seems costly
                assert(isPaired(next) && isPaired(s));  // double check?

                stack<implicitS> st_pairs;
                st_pairs.push(next);
                int cur_v = cv;
                while (!st_pairs.empty()) {
                    implicitS top = st_pairs.top();
                    st_pairs.pop();
                    assert(top.getConstVertices().size() == 2);  // ensure it is an edge (saddle)
                    cur_v = top.getConstVertices()[0] == cur_v ? top.getConstVertices()[1] : top.getConstVertices()[0];
                    vpath.push_back(cur_v);
                    implicitS ns(cur_v);
                    implicitS nnext;
                    if (getPair(ns, nnext)) {
                        assert(isPaired(nnext) && isPaired(ns));                 // double check
                        assert(nnext.getDim() == 1 && !sc.theSame(nnext, top));  // 1 is from stCell.getDim()
                        st_pairs.push(nnext);
                    }
                    //                    else {
                    //                        vids.push_back(ns.getConstVertices()[0]);
                    //                    }
                }
            }
            //            else {// critical 0-cell is on the 1-cell edge
            //                vids.push_back(s.getConstVertices()[0]);
            //            }

            v_paths.push_back(vpath);
        }
    }

    // organize cells
    // consider fork cases... -> one saddles has two same critical vertices, some
    // sep vpaths have some v-paths in common
    map<int, int> old2new;  // original point id -> written vtk point id
    int pts_num = 0;
    size_t cell_num = v_paths.size();
    size_t cell_items_num = 0;
    for (const auto &vp: v_paths) {
        cell_items_num = cell_items_num + 1 + vp.size();
        for (auto &v: vp) {
            if (old2new.find(v) == old2new.end()) {
                // new point in old2new
                old2new[v] = pts_num;
                pts_num += 1;
            }
        }
    }
    vector<int> new2old(static_cast<unsigned long>(pts_num));
    for (const auto &v: old2new) {
        new2old[v.second] = v.first;
    }
    // write vtk file
    cout << "ready to write the vtk file: " << vtkfile << endl;
    cout << "x,y extent: " << sc.sc_min_x << ", " << sc.sc_min_y << endl;

    ofstream ofs(vtkfile);
    ofs << fixed << setprecision(2);
    ofs << "# vtk DataFile Version 2.0\n"
           "tm graph\n"
           "ASCII\n"
           "DATASET UNSTRUCTURED_GRID\n";
    // points
    ofs << "POINTS " << pts_num << " float\n";
    for (const auto &v: new2old) {
        Vertex ver = sc.getVertex(v);
        if (scaled) {
            ofs << ver.getCoordinate(0) - sc.sc_min_x << " " << ver.getCoordinate(1) - sc.sc_min_y << " "
                << ver.getCoordinate(2) << endl;
        } else {
            ofs << ver.getCoordinate(0) << " " << ver.getCoordinate(1) << " " << ver.getCoordinate(2) << endl;
        }
    }
    // cells: polylines = 4
    ofs << "CELLS " << cell_num << " " << cell_items_num << endl;
    for (const auto &vp: v_paths) {
        ofs << vp.size() << " ";
        for (auto &v: vp) {
            ofs << old2new[v] << " ";
        }
        ofs << endl;
    }

    ofs << "CELL_TYPES " << cell_num << endl;
    for (int i = 0; i < cell_num; ++i) {
        ofs << "4 ";
    }
    ofs << endl;

    // field data
    ofs << "CELL_DATA " << cell_num << endl;
    ofs << "FIELD FieldData 1\n";
    ofs << "pid 1 " << cell_num << " int\n";
    for (int i = 0; i < cell_num; ++i) {
        ofs << i << " ";
    }
    ofs << endl;
    ofs.close();

    // debug
    if (outdebugfile) {
        ofstream df("debug.txt");
        df << "sep vpaths\n";
        size_t sp_id = 0;
        for (const auto &vp: v_paths) {
            df << sp_id << ": ";
            for (auto &v: vp) {
                df << v << " ";
            }
            sp_id += 1;
            df << endl;
        }
        df << "ol2new\n";
        for (auto on: old2new) {
            df << on.first << " <-> " << on.second << endl;
        }
        df.close();
    }
}

void FormanGradient::xx_vis_VPath_12(const string &vtkfile, const bool &scaled, const bool &outdebugfile) {

    std::vector<std::vector<implicitS>> vpaths;
    for (const auto &cc_1: criticalS[1]) {  // critical 1-cells
        //cout<<"\ncheck saddle: "<< cc_1<<endl;
        std::vector<std::vector<implicitS>> tmp_ps;
        saddle2max_vpaths(cc_1, tmp_ps);
        vpaths.insert(vpaths.end(), tmp_ps.begin(), tmp_ps.end());
    }

    std::cout << "\nvpaths#: " << vpaths.size() << endl;

    // todo: write to vtk
    ofstream ofs(vtkfile);
    ofs << fixed << setprecision(2);
    ofs << "# vtk DataFile Version 2.0\n"
           "tm graph\n"
           "ASCII\n"
           "DATASET UNSTRUCTURED_GRID\n";

    size_t pts_num = 0;
    size_t cell_num = 0;
    size_t cell_items_num = 0;
    for (const auto &vp: vpaths) {
        cell_num += vp.size();
        for (const auto &s: vp) {
            cell_items_num = cell_items_num + (1 + s.getConstVertices().size());
            pts_num += s.getConstVertices().size();
        }
    }

    ofs << "POINTS " << pts_num << " float\n";
    for (const auto &vp: vpaths) {
        for (const auto &s: vp) {
            for (const auto &vid: s.getConstVertices()) {
                Vertex ver = sc.getVertex(vid);
                if (scaled) {
                    ofs << ver.getCoordinate(0) - sc.sc_min_x << " " << ver.getCoordinate(1) - sc.sc_min_y << " "
                        << ver.getCoordinate(2) << endl;
                } else {
                    ofs << ver.getCoordinate(0) << " " << ver.getCoordinate(1) << " " << ver.getCoordinate(2) << endl;
                }
            }
        }
    }

    ofs << "CELLS " << cell_num << " " << cell_items_num << endl;
    int vid = 0;
    for (const auto &vp: vpaths) {
        for (const auto &s: vp) {
            size_t tmp_vnum = s.getConstVertices().size();
            if (tmp_vnum == 2) {
                ofs << "2 "; // VTK_LINE (=3)
                ofs << vid << " " << vid + 1 << endl;
                vid = vid + 2;
            }
            if (tmp_vnum == 3) {
                ofs << "3 "; //VTK_TRIANGLE(=5)
                ofs << vid << " " << vid + 1 << " " << vid + 2 << endl;
                vid = vid + 3;
            }
        }
    }

    ofs << "CELL_TYPES " << cell_num << endl;
    for (const auto &vp: vpaths) {
        for (const auto &s: vp) {
            size_t tmp_vnum = s.getConstVertices().size();
            if (tmp_vnum == 2) {
                ofs << "3 "; // VTK_LINE (=3)
            }
            if (tmp_vnum == 3) {
                ofs << "5 "; //VTK_TRIANGLE(=5)
            }
        }
    }
    ofs << endl;

    // field data
    ofs << "CELL_DATA " << cell_num << endl;
    ofs << "FIELD FieldData 1\n";
    ofs << "pid 1 " << cell_num << " int\n";
    int vp_id = 0;
    for (const auto &vp: vpaths) {
        for (const auto &s: vp) {
            ofs << vp_id << " ";
        }
        vp_id = vp_id + 1;
    }
    ofs << endl;
    ofs.close();

    ofs.close();

}

void FormanGradient::xx_visMorse(const string &output_prefix) {
    list<SSet> descending1manifold;
    list<SSet> ascending1manifold;
    for (auto criticalLVL: criticalS) {
        for (implicitS c: criticalLVL.second) {
            SSet ascending, descending;

            if (c.getDim() == 1) {
                computeDescendingCell(true, c, descending);
                descending1manifold.push_back(descending);

                computeAscendingCell(true, c, ascending);
                ascending1manifold.push_back(ascending);

            }
        }
    }

    string outfile1, outfile2;
    outfile1 = output_prefix + "_descending1cells.vtk";
    outfile2 = output_prefix + "_ascending1cells.vtk";
    print_out(outfile1.c_str(), descending1manifold, 3, 2);
    accurate_asc1cells(ascending1manifold, outfile2);

}

void FormanGradient::xx_critical_cells_net_2d(const std::string &vtkfile, const std::string &txtfile) {
    std::map<implicitS, int> cps_map;
    int cp_id = 0;
    std::vector<std::vector<int>> edges;
    for (const auto &s: criticalS[1]) {
        cps_map[s] = cp_id;
        cp_id = cp_id + 1;
        vector<implicitS> mins, maxs;
        saddle2minima(s, mins);
        saddle2maxima(s, maxs);
        for (const auto &m: mins) {
            if (cps_map.find(m) == cps_map.end()) {
                cps_map[m] = cp_id;
                cp_id = cp_id + 1;
            }
            int v1, v2;
            v1 = cps_map[s];
            v2 = cps_map[m];
            edges.push_back({v1, v2});
        }
        for (const auto &m: maxs) {
            if (cps_map.find(m) == cps_map.end()) {
                cps_map[m] = cp_id;
                cp_id = cp_id + 1;
            }
            int v1, v2;
            v1 = cps_map[s];
            v2 = cps_map[m];
            edges.push_back({v1, v2});
        }
    }

    std::vector<implicitS> cps(cps_map.size());
    for (const auto &c: cps_map) {
        cps[c.second] = c.first;
    }

    size_t pts_num, cell_num;
    pts_num = cps.size();
    cell_num = edges.size();

    cout << "pts#: " << pts_num << ", cells#: " << cell_num << endl;

    // write vtk
    ofstream ofs(vtkfile);
    ofs << fixed << setprecision(3);
    ofs << "# vtk DataFile Version 2.0\n"
           "critical cells\n"
           "ASCII\n"
           "DATASET UNSTRUCTURED_GRID\n";
    // points
    ofs << "POINTS " << pts_num << " float\n";
    for (const auto &c: cps) {
        Vertex v = sc.barycenter(c);
        ofs << v.getCoordinate(0) << " " << v.getCoordinate(1) << " " << v.getCoordinate(2) << endl;
    }


    // cells: VTK_LINE (=3)
    ofs << "CELLS " << cell_num << " " << cell_num * 3 << endl;
    for (const auto &e: edges) {
        ofs << "2 " << e[0] << " " << e[1] << endl;
    }

    ofs << "CELL_TYPES " << cell_num << endl;
    for (int i = 0; i < cell_num; ++i) {
        ofs << "3 "; //VTK_LINE (=3)
    }
    ofs << endl;

    // field data
    ofs << "CELL_DATA " << cell_num << endl;
    ofs << "FIELD FieldData 1\n";
    ofs << "edge_type 1 " << cell_num << " int\n";
    for (const auto &e: edges) {
        int sid1, sid2; // critical cell id, based on cps vector
        sid1 = e[0];
        sid2 = e[1];
        int d1, d2;
        d1 = cps[sid1].getDim();
        d2 = cps[sid2].getDim();
        if (d1 == 0 || d2 == 0) {
            ofs << "0 ";
        }
        if (d1 == 2 || d2 == 2) {
            ofs << "2 ";
        }
    }
    ofs << endl;
    ofs.close();


    //  write txt file
    ofstream ofs1(txtfile);
    ofs1 << fixed << setprecision(3);
    ofs1 << pts_num << " " << cell_num << endl;
    // points: x y z dim
    for (const auto &c: cps) {
        Vertex v = sc.barycenter(c);
        ofs1 << v.getCoordinate(0) << " " << v.getCoordinate(1) << " " << v.getCoordinate(2) << " " << c.getDim()
             << endl;

    }
    // edges
    for (const auto &e: edges) {
        ofs1 << e[0] << " " << e[1] << endl;
    }
    ofs1.close();

}

void FormanGradient::writeVTK_gradient_2d(const string &vtkfile) {
    // todo: check
    // https://github.com/IuricichF/FormanGradient2D/blob/master/source/LibForman/io.cpp#L221

    // point to find its paired edge
    // edge to find its paired triangle or use triangle to find paired edge,
    // as point and triangle are explictly encoded in IA?

    cout << "\noutput gradient 2d\n";

    std::vector<std::array<float, 3>> v_coords, v_gradients;
    for (auto d: sc.getTopSimplexesSet()) {
        cout << "dim: " << d << endl;
    }
    // check vertex -> edge
    for (int ts_vid = 0; ts_vid < sc.getVerticesNum(); ++ts_vid) {
        //for (auto ts: sc.getTopSimplices(0)) {
        //int ts_vid = ts.getVertexIndex(0);
        vector<float> ts_coords = sc.getVertex(ts_vid).getCoordinates();
        v_coords.push_back({ts_coords[0], ts_coords[1], ts_coords[2]});

        implicitS s(ts_vid);  // vertex simplex
        implicitS next;
        // find paired edge
        if (getPair(s, next)) {
            // need?
            // assert(isPaired(next) && isPaired(s));  // double check
            // get the other point on point
            vector<int> next_vids = next.getConstVertices();
            int tmp_vid = next_vids[0];
            if (tmp_vid == ts_vid) {
                tmp_vid = next_vids[1];
            }
            vector<float> tmp_coords = sc.getVertex(tmp_vid).getCoordinates();
            float del_x, del_y, del_z;
            del_x = tmp_coords[0] - ts_coords[0];
            del_y = tmp_coords[1] - ts_coords[1];
            del_z = tmp_coords[2] - ts_coords[2];
            v_gradients.push_back({del_x, del_y, del_z});
        } else { v_gradients.push_back({0.0, 0.0, 0.0}); }
    }

    // todo: edge to triangles
    //  check the code.

    // debug write vtk
    size_t pts_num = v_gradients.size();
    ofstream ofs(vtkfile);
    ofs << fixed << setprecision(3);
    ofs << "# vtk DataFile Version 2.0\n"
           "critical cells\n"
           "ASCII\n"
           "DATASET UNSTRUCTURED_GRID\n";
    // points
    ofs << "POINTS " << pts_num << " float\n";
    for (const auto &pt: v_coords) {
        ofs << pt[0] << " " << pt[1] << " " << pt[2] << endl;
    }

    ofs << "POINT_DATA " << pts_num << endl;
    ofs << "VECTORS vector float\n";
    for (const auto &g: v_gradients) {
        ofs << g[0] << " " << g[1] << " " << g[2] << endl;
    }
    ofs.close();
}

void FormanGradient::write_critical_persistence_2d(const std::string &outfile) {
    // NOT ready to use.
    cout << "\nwrite critical persistence infos\n";
    std::map<implicitS, std::array<float, 2>> min_pers, saddle_pers, max_pers;
    // mins
    for (const auto &m: criticalS[0]) {
        float val = sc.getVertex(m.getConstVertices()[0]).getCoordinate(2);
        min_pers[m][0] = val;

    }
    // saddle to mins
    for (const auto &saddle: criticalS[1]) {  // edge saddles
        float tmp_val1, tmp_val2, val;
        tmp_val1 = sc.getVertex(saddle.getConstVertices()[0]).getCoordinate(2);
        tmp_val2 = sc.getVertex(saddle.getConstVertices()[1]).getCoordinate(2);
        val = max(tmp_val1, tmp_val2);
        vector<implicitS> minima;
        saddle2minima(saddle, minima);
        implicitS m1, m2;
        m1 = minima[0];
        m2 = minima[1];
        min_pers[m1][1] = val;
        min_pers[m2][1] = val;
        saddle_pers[saddle][0] = val;

    }

    // note: 1 saddle -> 2 mins, 2 maxs
    std::map<implicitS, int> s2nums;
    // max to saddle
    for (const auto &m: criticalS[2]) {  //
        float tmp_val1, tmp_val2, tmp_val3, val;
        tmp_val1 = sc.getVertex(m.getConstVertices()[0]).getCoordinate(2);
        tmp_val2 = sc.getVertex(m.getConstVertices()[1]).getCoordinate(2);
        tmp_val3 = sc.getVertex(m.getConstVertices()[2]).getCoordinate(2);
        val = max(tmp_val1, tmp_val2);
        val = max(val, tmp_val3);
        vector<implicitS> saddles;
        // todo: finish it
        //max2saddle(m, saddles);
        // debug: 1 max -> varied number of saddles
        // cout << "matched critical saddles#: " << saddles.size() << endl;
        //
        for (const auto &s: saddles) {
            saddle_pers[s][1] = val;
            s2nums[s]++;
        }
        max_pers[m][0] = val;
        max_pers[m][1] = sc.sc_max_z;
    }


    // todo: debug: check 1 saddle match to how many max. some saddle only found 1 max
    cout << "check saddle -> max#. should be 2\n";
    for (const auto &s: s2nums) {
        implicitS tmp_s = s.first;
        int tmp_num = s.second;
        if (s.second != 2) {
            cout << tmp_s << ": " << tmp_num << endl;
        }
    }
    size_t s_num = saddle_pers.size();


    // write output
    // x1 y1 z1 x2 y2 z2 ... v1 v2
    ofstream ofs(outfile);
    ofs << fixed << setprecision(3);
    ofs << min_pers.size() << " " << s_num << " " << max_pers.size() << endl;
    for (const auto &k: min_pers) {
        implicitS tmp_s = k.first;
        std::vector<float> coords = sc.getVertex(tmp_s.getConstVertices()[0]).getCoordinates();
        std::array<float, 2> tmp = k.second;
        ofs << coords[0] << " " << coords[1] << " " << coords[2] << " " << tmp[0] << " " << tmp[1] << endl;
    }
    for (const auto &k: saddle_pers) {
        implicitS tmp_s = k.first;
        std::vector<float> coords1 = sc.getVertex(tmp_s.getConstVertices()[0]).getCoordinates();
        std::vector<float> coords2 = sc.getVertex(tmp_s.getConstVertices()[1]).getCoordinates();
        std::array<float, 2> tmp = k.second;
        ofs << coords1[0] << " " << coords1[1] << " " << coords1[2] << " ";
        ofs << coords2[0] << " " << coords2[1] << " " << coords2[2] << " ";
        ofs << tmp[0] << " " << tmp[1] << endl;
    }
    for (const auto &k: max_pers) {
        implicitS tmp_s = k.first;
        std::vector<float> coords1 = sc.getVertex(tmp_s.getConstVertices()[0]).getCoordinates();
        std::vector<float> coords2 = sc.getVertex(tmp_s.getConstVertices()[1]).getCoordinates();
        std::vector<float> coords3 = sc.getVertex(tmp_s.getConstVertices()[2]).getCoordinates();

        ofs << coords1[0] << " " << coords1[1] << " " << coords1[2] << " ";
        ofs << coords2[0] << " " << coords2[1] << " " << coords2[2] << " ";
        ofs << coords3[0] << " " << coords3[1] << " " << coords3[2] << " ";
        std::array<float, 2> tmp = k.second;
        ofs << tmp[0] << " " << tmp[1] << endl;
    }
    ofs.close();
}