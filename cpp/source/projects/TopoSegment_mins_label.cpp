#include "TopoSegment.h"


int TopoSegment::initSeeds_mins_by_as() {
    // init _min2lbl: min: lbl (-1,-2, ...)
    if (_showlog)
        cout << "init _min2lbl\n";
    std::vector<implicitS> all_mins(criticalS[0].begin(), criticalS[0].end()); // all mins.
    for (int i = 0; i < all_mins.size(); ++i) {
        _min2lbl[all_mins[i]] = -1 - i;
    }

    // find connected points in the low height range,
    if (_showlog)
        cout << "find connected mins in the low height\n";
    std::vector<std::set<int> > con_low_pts;
    // Note: it iterates all points, try to improve the performance
    _get_connected_points(sc.sc_min_z, sc.sc_min_z + 0.5, con_low_pts);
    if (_showlog) {
        cout << "# of connected objects: " << con_low_pts.size() << endl;
        for (const auto& g : con_low_pts) {
            cout << "group: ";
            cout << "#: " << g.size();
            cout << endl;
        }
        cout << endl;
    }

    // create a map to improve the speed
    // map: low point id: label
    std::map<int, int> lowpt2lbl; // lbl starts from 1
    for (int i = 0; i < con_low_pts.size(); ++i) {
        std::set<int> tmp = con_low_pts[i];
        for (const auto& pid : tmp) {
            lowpt2lbl[pid] = i + 1; // lbl starts from 1
        }
    }

    // low mins are seeds. connected mins have the same label
    for (const auto& m : all_mins) {
        int vid = m.getConstVertices()[0];
        Vertex v = sc.getVertex(vid);
        if (v.getCoordinate(2) > sc.sc_min_z + 0.5) {
            // xyz
            continue;
        }
        if (_showlog)
            cout << "check min: " << vid << endl;
        int tmp_lbl = lowpt2lbl[vid];
        _min2lbl[m] = tmp_lbl;
        if (_showlog)
            cout << "-> set as seed with id: " << tmp_lbl << endl;
    } // build seed_mins
    if (_showlog) {
        cout << "# of stems: " << con_low_pts.size() << endl;
    }

    return 0;
}


int TopoSegment::initSeeds_mins_by_trunkpts_xyz(
    const std::vector<std::vector<float> >& seeds_pts,
    const double& th_p2trunk_distance,
    const double& th_search_radius) {
    /* key ideas:
     * default: th_p2trunk_distance = 0.2,  th_search_radius = 0.25
     * find initial labeled mins: close to seeds.
     *  distance: min to seed pts or min to seed centroid
     *  min.z between trunk info zmin, zmax
     *  distance (min, trunk_pts) < th_p2trunk_distance
     * if no exist, then
     *  distance (min, trunk_cp_pt) < th_search_radius
     *
     * note: valid label starts from 1.
     */

    // init _min2lbl: min: lbl (-1,-2, ...)
    if (_showlog)
        cout << "\ninit _min2lbl\n";
    std::vector<implicitS> all_mins(criticalS[0].begin(), criticalS[0].end()); // all mins.
    for (int i = 0; i < all_mins.size(); ++i) {
        _min2lbl[all_mins[i]] = -1 - i;
    }

    // // read seed pts (x,y, ...) from seed_file
    // std::vector<std::vector<float> > seeds_pts;
    // tts_read_pts_fs(seed_file, seeds_pts);

    std::map<int, std::set<std::array<float, 3> > > trunk_pts; // todo: use seed2pids.
    for (const auto& pt : seeds_pts) {
        std::array<float, 3> tmp_xyz = {pt[0], pt[1], pt[2]};
        int tmp_l = (int)pt[3];
        trunk_pts[tmp_l].insert(tmp_xyz);
    }
    if (trunk_pts.empty()) {
        cout << "no trunk provided. use alpha shape components as default\n";
        return 1; // something wrong
    }
    if (_showlog)
        std::cout << "seeds #: " << trunk_pts.size() << std::endl;

    // compute trunk infos: x,y,zmin,zmax
    std::map<int, std::array<float, 4> > trunk_infos; // (x,y,zmin,zmax) for each trunk
    for (const auto& tt : trunk_pts) {
        //std::cout << "input seed label: " << tt.first << ", seed#: " << tt.second.size() << endl;
        int trunk_id = tt.first;
        std::set<std::array<float, 3> > tt_xyz = tt.second;
        float cx, cy, tmp_zmax, tmp_zmin;
        cx = 0.0;
        cy = 0.0;
        tmp_zmax = std::numeric_limits<float>::min();
        tmp_zmin = std::numeric_limits<float>::max();
        for (const auto& tmp_p : tt_xyz) {
            //std::cout << tmp_p[0] << " " << tmp_p[1] << endl;
            cx = cx + tmp_p[0];
            cy = cy + tmp_p[1];
            tmp_zmin = min(tmp_p[2], tmp_zmin);
            tmp_zmax = max(tmp_p[2], tmp_zmax);
        }
        //std::cout << "sum" << ": " << "(x,y): " << cx << ", " << cy << endl;
        //std::cout << "size:" << tt_xyz.size() << endl;
        cx = cx / (tt_xyz.size());
        cy = cy / (tt_xyz.size());
        //std::cout << tt.first << ": " << "(x,y): " << cx << ", " << cy << endl;
        trunk_infos[trunk_id] = {cx, cy, tmp_zmin, tmp_zmax};
    }

    // find initial labeled mins: close to seeds.
    // distance: min to seed pts or min to seed centroid
    // min.z between trunk info zmin, zmax
    // distance (min, trunk_pts) should be smaller than th_p2trunk_distance

    for (const auto& m : all_mins) {
        int vid = m.getConstVertices()[0];
        Vertex v = sc.getVertex(vid);
        // here, check mins from all heights ranges
        // note: we now assume trunks are straight, we may have a better way to define the distance from min to a trunk
        if (_showlog)
            cout << "check min: " << vid << endl;
        std::vector<float> v_coords = v.getCoordinates();
        std::array<float, 3> v_xyz = {v_coords[0], v_coords[1], v_coords[2]};

        // find the closet trunk to min. check distance to the trunk/loc
        // distance: p to the trunk pts (not the center pts)
        int bst_trunk_lbl = 0;
        float tmp_min_dis_sq = std::numeric_limits<float>::max(); //init dis. large
        for (const auto& tmp_tloc : trunk_pts) {
            // check trunk first
            // valid seed/trunk: length >= 2.0m.
            // min.z need to be [trunk.zmin, trunk.zmax]
            int tmp_trunk_id = tmp_tloc.first;
            std::array<float, 4> t_infos = trunk_infos[tmp_trunk_id]; // (cx, cy, zmin, zmax)
            // no need to consider the trunk length now.
            // if ((t_infos[3] - t_infos[2]) < 2.0) { continue; }
            // min.z need to be [trunk.zmin, trunk.zmax]
            if ((v_coords[2] > t_infos[3]) || (v_coords[2] < t_infos[2])) {
                continue;
            }
            std::set<std::array<float, 3> > tmp_trunk_pts = tmp_tloc.second;
            // min to trunk pts minimal distance
            float tmp_d_sq = dis_sq_pt2pts(v_xyz, tmp_trunk_pts);
            if (tmp_d_sq < th_p2trunk_distance * th_p2trunk_distance) {
                if (tmp_d_sq < tmp_min_dis_sq) {
                    tmp_min_dis_sq = tmp_d_sq;
                    bst_trunk_lbl = tmp_trunk_id;
                }
            }
        } // find the closest trunk
        if (bst_trunk_lbl == 0) {
            // not found the closest trunk
            // to find the closest trunk bottoms (distance: centroid point)
            // note here we work on xyz
            for (const auto& tmp_tloc : trunk_pts) {
                int tmp_trunk_id = tmp_tloc.first;
                std::array<float, 4> t_infos = trunk_infos[tmp_trunk_id]; // (cx,cy, zmin, zmax)
                // min.z need to be <= trunk.zmax. min should be < trunk.zmin
                if (v_coords[2] > t_infos[3]) { continue; }

                // min to centroid distance
                std::array<float, 2> t_coords = {t_infos[0], t_infos[1]}; // (x,y)
                // check 2d distance, in case trunk is not complete, we use trunk 2d locs to init nearby mins
                float tmp_d2_sq = (v_coords[0] - t_coords[0]) * (v_coords[0] - t_coords[0]) +
                                  (v_coords[1] - t_coords[1]) * (v_coords[1] - t_coords[1]);
                if (tmp_d2_sq < th_search_radius * th_search_radius) {
                    if (tmp_d2_sq < tmp_min_dis_sq) {
                        tmp_min_dis_sq = tmp_d2_sq;
                        bst_trunk_lbl = tmp_trunk_id;
                    }
                }
            }
        } // find the closest bottom

        // Note: don't assign 0 to the min. 0 means unlabeled, but low min should be labeled.
        if (bst_trunk_lbl != 0) {
            _min2lbl[m] = bst_trunk_lbl;
            if (_showlog)
                cout << "-> set as seed with id: " << bst_trunk_lbl << endl;
        }
    } // build seed_mins

    // debug
    if (_showlog) {
        map<int, std::vector<implicitS> > gp_mins;
        for (const auto& ml : _min2lbl) {
            gp_mins[ml.second].push_back(ml.first);
        }
        // output mins with labels
        _output_gpmins_pts(gp_mins, _workspace + "/" + _file_name + "_test_mins_seed.pts", false);
    }

    return 0;
}


int TopoSegment::growFromSeeds_mins_basic() {
    // low mins are seeds. connected mins have the same label
    std::vector<implicitS> sorted_seed_mins; // low 2 high z
    // todo: [long term] make sorted_seed_mins like automated sorted vector?
    for (const auto& m : _min2lbl) {
        if (m.second > 0) {
            sorted_seed_mins.push_back(m.first);
        }
    }
    // sort seeds
    std::sort(sorted_seed_mins.begin(), sorted_seed_mins.end(), [&](const implicitS& m1, const implicitS& m2) {
        Vertex v1, v2;
        v1 = sc.getVertex(m1.getConstVertices()[0]);
        v2 = sc.getVertex(m2.getConstVertices()[0]);
        return v1.getCoordinate(2) < v2.getCoordinate(2); // xyz
    });
    if (_showlog) { cout << "# of seeds: " << sorted_seed_mins.size() << endl; }

    // all mins are seeds, no need to merge
    std::vector<implicitS> all_mins(criticalS[0].begin(), criticalS[0].end()); // all mins.
    if ((all_mins.size() - sorted_seed_mins.size()) == 0) {
        if (_showlog)
            cout << "all mins are seeds.\n";
        // all mins are seed mins, no isolated mins, no non-seed mins
        return 0;
    } // all mins are seeds, no need to merge

    // create min2mins
    std::map<implicitS, std::set<implicitS> > min2mins; // critical minima to its directed connected critical minima
    for (const auto& saddle : criticalS[1]) {
        // edge saddles
        vector<implicitS> minima;
        saddle2minima(saddle, minima);
        implicitS m1, m2;
        m1 = minima[0];
        m2 = minima[1];
        if (m1 == m2) {
            continue;
        } else {
            min2mins[m1].insert(m2);
            min2mins[m2].insert(m1);
        }
    } // create min2mins by saddles

    // label mins
    while (!sorted_seed_mins.empty()) {
        //  update the sorted seeds mins in each round
        if (_showlog)
            cout << "\nstart round\n";
        // sort seeds
        std::sort(sorted_seed_mins.begin(), sorted_seed_mins.end(), [&](const implicitS& m1, const implicitS& m2) {
            Vertex v1, v2;
            v1 = sc.getVertex(m1.getConstVertices()[0]);
            v2 = sc.getVertex(m2.getConstVertices()[0]);
            return v1.getCoordinate(2) < v2.getCoordinate(2); // xyz
        });

        // remove seed ids
        std::vector<int> rm_seed_ids;
        std::set<implicitS> add_seeds;
        // iterate seeds to grow
        // if one min grows, then try to label ALL its unlabeled neighbors
        for (int i = 0; i < sorted_seed_mins.size(); ++i) {
            // const auto& seed_min : sorted_seed_mins
            implicitS seed_min = sorted_seed_mins[i];
            rm_seed_ids.push_back(i);

            bool min_grows = false; // does this seed grow?
            int tmp_lbl = _min2lbl[seed_min];
            //std::array<float, 4> tmp_trunk_info = trunk_infos.at(tmp_lbl);
            for (const auto& nbh : min2mins[seed_min]) {
                int nbh_lbl = _min2lbl[nbh];
                if (nbh_lbl < 0) {
                    // unlabeled min
                    _min2lbl[nbh] = tmp_lbl;
                    add_seeds.insert(nbh);
                    min_grows = true;
                }
            }
            // if one seed is grown
            // stop checking other mins, this growing round is finished
            if (min_grows) {
                break;
            }
        } // iterate min in sorted_seed_mins

        // remove seeds from large id to small, for erase operation
        // reverse vector
        std::reverse(rm_seed_ids.begin(), rm_seed_ids.end());
        if (_showlog)
            cout << "remove processed seeds:\n";
        for (const auto& rm_id : rm_seed_ids) {
            if (_showlog)
                cout << rm_id << " ";
            sorted_seed_mins.erase(sorted_seed_mins.begin() + rm_id);
        }
        if (_showlog)
            cout << endl;

        // add seeds
        if (_showlog)
            cout << "add new seeds:\n";
        for (const auto& nbr : add_seeds) {
            if (_showlog)
                cout << nbr << " ";
            sorted_seed_mins.push_back(nbr);
        }
        if (_showlog)
            cout << endl;

        //
        if (_showlog) {
            cout << "after growing seeds:\n";
            cout << "# of seeds: " << sorted_seed_mins.size() << endl;
        }
    } // label rest labels

    return 0;
}


string TopoSegment::label_pts_from_mins(const bool& output) {
    if (_showlog)
        cout << "\nlabel points\n";
    // init
    _pts_lbls = std::vector<int>(sc.getVerticesNum(), 0);
    // label pts
    for (const auto& ml : _min2lbl) {
        implicitS cc = ml.first;
        int lbl = ml.second;
        SSet tmp_cells;
        computeAscendingCell(true, cc, tmp_cells);
        // label points in the ascending/descending cells
        for (const auto& ac : tmp_cells) {
            for (const auto& vid : ac.getConstVertices()) {
                if (_pts_lbls[vid] > 0) {
                    int old_lbl = _pts_lbls[vid];
                    if (old_lbl != lbl) {
                        _pts_lbls[vid] = 0; //  point is labeled at least twice
                    }
                } else {
                    _pts_lbls[vid] = lbl;
                }
            }
        }
    }

    // support text .off and binary .ply
    string outfile = _workspace + _file_name + "_lbl.pts";
    if (_infile.find(".off") != std::string::npos) {
        outfile = _workspace + _file_name + "_lbl.pts";
        if (output) {
            cout << "write: " << outfile << endl;
            _output_pts_with_label_pts(outfile, _pts_lbls, false);
        }
    } else if (_infile.find(".ply") != std::string::npos) {
        outfile = _workspace + _file_name + "_lbl.ply";
        if (output) {
            cout << "write: " << outfile << endl;
            _output_pts_with_label_pts_ply(outfile, _pts_lbls, false);
        }
    }

    return outfile;
}

string TopoSegment::label_pts_from_mins_parallel(const bool& output) {
    if (_showlog)
        cout << "\nlabel points\n";
    // init
    _pts_lbls = std::vector<int>(sc.getVerticesNum(), 0);

    // label pts
    // for (const auto &ml: _min2lbl) {
    //     implicitS cc = ml.first;
    //     int lbl = ml.second;
    //     SSet tmp_cells;
    //     computeAscendingCell(true, cc, tmp_cells);
    //     // label points in the ascending/descending cells
    //     for (const auto &ac: tmp_cells) {
    //         for (const auto &vid: ac.getConstVertices()) {
    //             if (_pts_lbls[vid] > 0) {
    //                 int old_lbl = _pts_lbls[vid];
    //                 if (old_lbl != lbl) {
    //                     _pts_lbls[vid] = 0; //  point is labeled at least twice
    //                 }
    //             } else {
    //                 _pts_lbls[vid] = lbl;
    //             }
    //         }
    //     }
    // }

    // Convert map to vector for better parallel iteration
    std::vector<std::pair<implicitS, int> > min2lbl_vec(_min2lbl.begin(), _min2lbl.end());

    // label pts
#pragma omp parallel num_threads(6)
    {
        // Thread-local vector for labels
        std::vector<int> local_labels(sc.getVerticesNum(), 0);

#pragma omp for schedule(dynamic) nowait
        for (size_t i = 0; i < min2lbl_vec.size(); ++i) {
            const auto& cc = min2lbl_vec[i].first;
            int lbl = min2lbl_vec[i].second;
            SSet tmp_cells;
            computeAscendingCell(true, cc, tmp_cells);

            // label points in the ascending/descending cells
            for (const auto& ac : tmp_cells) {
                for (const auto& vid : ac.getConstVertices()) {
                    if (local_labels[vid] > 0) {
                        int old_lbl = local_labels[vid];
                        if (old_lbl != lbl) {
                            local_labels[vid] = 0; //  point is labeled at least twice
                        }
                    } else {
                        local_labels[vid] = lbl;
                    }
                }
            }
        }

        // Merge results from all threads
#pragma omp critical
        {
            for (size_t i = 0; i < _pts_lbls.size(); ++i) {
                if (local_labels[i] != 0) {
                    if (_pts_lbls[i] > 0 && _pts_lbls[i] != local_labels[i]) {
                        _pts_lbls[i] = 0;
                    } else {
                        _pts_lbls[i] = local_labels[i];
                    }
                }
            }
        }
    }

    // support text .off and binary .ply
    string outfile = _workspace + _file_name + "_lbl.pts";
    if (_infile.find(".off") != std::string::npos) {
        outfile = _workspace + _file_name + "_lbl.pts";
        if (output) {
            cout << "write: " << outfile << endl;
            _output_pts_with_label_pts(outfile, _pts_lbls, false);
        }
    } else if (_infile.find(".ply") != std::string::npos) {
        outfile = _workspace + _file_name + "_lbl.ply";
        if (output) {
            cout << "write: " << outfile << endl;
            _output_pts_with_label_pts_ply(outfile, _pts_lbls, false);
        }
    }

    return outfile;
}

string TopoSegment::output_labeled_pts(const string& outfile) {
    string ret_str = outfile;
    if (outfile.empty()) {
        ret_str = _workspace + _file_name + "_lbl.pts";
        if (_infile.find(".off") != std::string::npos) {
            _output_pts_with_label_pts(ret_str, _pts_lbls, false);
        } else if (_infile.find(".ply") != std::string::npos) {
            ret_str = _workspace + _file_name + "_lbl.ply";
            _output_pts_with_label_pts_ply(ret_str, _pts_lbls, false);
        }
    } else {
        if (_infile.find(".off") != std::string::npos) {
            _output_pts_with_label_pts(ret_str, _pts_lbls, false);
        } else if (_infile.find(".ply") != std::string::npos) {
            _output_pts_with_label_pts_ply(ret_str, _pts_lbls, false);
        }
    }
    return ret_str;
}