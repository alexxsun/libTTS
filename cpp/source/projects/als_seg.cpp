#include "als_seg.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <queue>
#include <set>
#include <sstream>

using namespace std;

bool xx_debug = false;

int als_segment(const string& infile, double th_radius, double th_forest_pers, double th_height, bool th_height_funtype,
                double th_pers_H, double th_pers_I, const string& ds) {
    /*
     * Input:
     * infile: pts file
     * th_radius: radius of connecting close points. >0: manually set; otherwise auto-set
     * th_forest_pers: forest persistence value rate. [0,1]. default: 0.02
     * th_height: height range to connect point. rate: [0,1] >0: manually set, otherwise auto-set. default: 0.5
     * th_height_funtype: height function type. true: height function, false: inverse height function
     * th_pers_H, th_pers_I: persistence value rates for height and inverse height cases. [0,1]. >0 manually set,
     *  otherwise auto-set
     *
     * Note:
     * point label: valid labels >=0 (starts from 0). not labeled: -1, boundary: -2
     */

    vector<vector<double> > points;
    read_points(infile, points);
    cout << "read pts file: " << infile << endl;

    // double th_radius;
    vector<double> pts_density;
    vector<vector<int> > VV;
    connect_projected_pts(points, th_radius, pts_density, VV);
    cout << "connect projected points\n";

    // debug
    //  exit(1);

    th_forest_pers = th_forest_pers < 0 ? 0.02 : th_forest_pers;
    cout << "forest persistence rate: " << th_forest_pers << endl;
    vector<int> pts_lbls;
    map<int, vector<int> > tree_clusters; // cluster id: points ids
    segment_forest(th_forest_pers, pts_density, VV, pts_lbls, tree_clusters);
    cout << "segment forest into tree clusters\n";

    // auto-set parameters follows paper strictly
    if (th_height_funtype) {
        // stem focus.
        cout << "height function. stem focus\n";
        th_pers_H = th_pers_H < 0 ? 0.01 : th_pers_H;
        th_pers_I = th_pers_I < 0 ? 0.7 : th_pers_I;
    } else {
        // treetop focus.
        cout << "inverse height function. treetop focus\n";
        th_pers_I = th_pers_I < 0 ? 0.01 : th_pers_I;
        th_pers_H = th_pers_H < 0 ? 0.7 : th_pers_H;
    }
    cout << "persistence of height and inverse height functions: " << th_pers_H << ", " << th_pers_I << endl;

    int regionsH, regionsI; // tree label. starts from 0
    regionsH = regionsI = 0;
    vector<int> pts_lbls_H, pts_lbls_I;
    pts_lbls_H = vector<int>(points.size(), -1); // init segment, which is the return value
    pts_lbls_I = vector<int>(points.size(), -1); // init segment, which is the return value
    th_height = th_height < 0 ? 0.5 : th_height;
    for (const auto& tc : tree_clusters) {
        vector<int> tc_pids = tc.second;
        vector<vector<int> > tc_VV;
        double zmin, zmax;
        compute_cylinder(th_height, th_radius, points, tc_pids, tc_VV, zmin, zmax);
        segment_cluster(points, zmin, zmax, tc_pids, tc_VV, true, th_pers_H, pts_lbls_H, regionsH);
        segment_cluster(points, zmin, zmax, tc_pids, tc_VV, false, th_pers_I, pts_lbls_I, regionsI);
    }
    cout << "segment tree clusters from top-down and bottom-up\n";

    // write each type of segmentation results
    string outfile1 = infile.substr(0, infile.size() - 4) + "_H_lbls.pts";
    output_labeled_pts(points, pts_lbls_H, outfile1);
    cout << "write: " << outfile1 << endl;
    outfile1 = infile.substr(0, infile.size() - 4) + "_I_lbls.pts";
    output_labeled_pts(points, pts_lbls_I, outfile1);
    cout << "write: " << outfile1 << endl;

    pts_lbls.clear();
    int trees_num = 0;
    merge_labels(pts_lbls_H, pts_lbls_I, pts_lbls, trees_num);
    cout << "merge segmentation results\n";
    if (ds.size() < 1 || ds == "none") {
        cout << "no post-process\n";
    } else {
        post_process(pts_lbls, trees_num, ds);
        cout << "post-process results\n";
    }

    string outfile = infile.substr(0, infile.size() - 4) + "_lbls.pts";
    output_labeled_pts(points, pts_lbls, outfile);
    cout << "write: " << outfile << endl;

    return 0;
}

template <typename T>
void xx_str2vector(const string& str, vector<T>& vec) {
    istringstream iss(str);
    // https://stackoverflow.com/questions/9986091/how-do-you-convert-a-string-into-an-array-of-doubles
    std::copy(std::istream_iterator<T>(iss), std::istream_iterator<T>(), std::back_inserter(vec));
}

int read_points(const string& ptsfile, vector<vector<double> >& pts) {
    ifstream ifs(ptsfile);
    string line;
    while (getline(ifs, line)) {
        std::vector<double> coord;
        xx_str2vector<double>(line, coord);
        pts.push_back(coord);
    }
    ifs.close();
    return 1;
}

int connect_projected_pts(const vector<vector<double> >& pts, double& Radius, vector<double>& density,
                          vector<vector<int> >& VV) {
    /* Inputs:
     * pts: all points
     * Outputs:
     * Radius
     * density: point density
     * VV: connected pids
     */
    // generate a graph of points projected on the ground
    // todo: delete ann pointers?
    typedef ANNcoord* ANNpoint;
    typedef ANNpoint* ANNpointArray;
    typedef ANNdist* ANNdistArray;

    typedef ANNidx* ANNidxArray;

    int vertexNumber = pts.size();
    ANNpointArray dataPts = annAllocPts(vertexNumber, 2);
    for (int i = 0; i < vertexNumber; i++) {
        ANNpoint newPt = annAllocPt(2);
        newPt[0] = pts[i][0];
        newPt[1] = pts[i][1];
        dataPts[i] = newPt;
    }

    ANNkd_tree* kdTree = new ANNkd_tree(dataPts, vertexNumber, 2);

    // search a good radius for the search
    if (Radius <= 0) {
        // auto-set Radius
        /*
         * old way: starts from 0.25, avgVV (50,150) ok, < 50: +0.25, otherwise -0.25
         * new way: starts from 0.5, max 1.0, avgVV<70, +0.1
         */
        Radius = 0.5; // 0.25;
        double avgVV = 0;
        bool vvwrong = true;
        while (vvwrong) {
            if (Radius >= 0.99) {
                break;
            }
            ANNpoint queryPt = annAllocPt(2);
            // start here
            for (int i = 0; i < vertexNumber; i++) {
                // double x, y, z;
                queryPt[0] = pts[i][0];
                queryPt[1] = pts[i][1];
                ANNidxArray ids = new ANNidx[0];
                ANNdistArray dists = new ANNdist[0];
                ANNdist dist = (double)pow(Radius, 2);
                int knn = kdTree->annkFRSearch(queryPt, dist, 0, ids, dists, 0.0);
                avgVV += (double)knn;
            }
            avgVV = avgVV / (double)vertexNumber;
            if (avgVV < 70) {
                Radius = Radius + 0.1;
            } else {
                vvwrong = false;
            }
            // old way
            //      if (avgVV < 150.0 && avgVV > 50.0) {
            //        vvwrong = false;
            //      } else if (avgVV < 50.0) {
            //        Radius = Radius + 0.25;
            //      } else {
            //        Radius = Radius - 0.25;
            //      }
        }
    } // auto-set Radius
    cout << "Selected Radius: " << Radius << endl;

    // build up the structure -- the graph of 2d points (projected on the xy plane)
    density = vector<double>(vertexNumber, 0); // density is the height value
    VV = vector<vector<int> >(vertexNumber, vector<int>());
    ANNpoint queryPt = annAllocPt(2);
    for (int i = 0; i < vertexNumber; i++) {
        queryPt[0] = pts[i][0];
        queryPt[1] = pts[i][1];
        density[i] = pts[i][2]; // init density as height z
        ANNidxArray ids = new ANNidx[0];
        ANNdistArray dists = new ANNdist[0];
        ANNdist dist = (double)pow(Radius, 2);
        int knn = kdTree->annkFRSearch(queryPt, dist, 0, ids, dists, 0.0);
        ids = new ANNidx[knn];
        dists = new ANNdist[knn];
        knn = kdTree->annkFRSearch(queryPt, dist, knn, ids, dists, 0.0);
        for (int j = 0; j < knn; j++) {
            if (ids[j] > i) {
                VV[i].push_back(ids[j]);
                VV[ids[j]].push_back(i);
            }
        }
    } // built VV and calculate density values

    // debug: find average number of neighbor
    if (Radius > 0) {
        double avg_nbr_num = 0;
        for (const auto& nbrs : VV) {
            avg_nbr_num = avg_nbr_num + nbrs.size();
        }
        avg_nbr_num = avg_nbr_num / VV.size();
        cout << "Radius: " << Radius << ", average nbrs #: " << avg_nbr_num << endl;
    }

    // 10 times of Laplacian_smoothing
    for (int i = 0; i < 10; i++) {
        vector<double> newDensity(vertexNumber, 0);
        for (int j = 0; j < vertexNumber; j++) {
            for (auto v : VV[j]) {
                newDensity[j] += density[v];
            }
            newDensity[j] /= (double)VV[j].size();
        }
        density = newDensity;
    }

    // reverse the density (height) values
    for (int i = 0; i < vertexNumber; i++) {
        density[i] = -density[i];
    }
    return 1;
}

int segment_forest(const double& perc_pers, const vector<double>& density, vector<vector<int> >& VV,
                   vector<int>& pts_lbls, map<int, vector<int> >& tree_clusters) {
    /*
     * given the graph of 2D points, segment the graph based on the watershed.
     * Inputs:
     * perc_pers: (0-1), real persistence value = std::abs(maxD - minD) * perc_pers
     * density: values of vertices
     * VV: the graph of vertices
     * Outputs:
     * pts_lbls: label of each vertex.
     * tree cluster:
     *
     */

    int totVertices = density.size();
    map<int, int> regionToMinima; // region id to root vertex id
    vector<int> labels(totVertices, -1); //  unlabeled: -1; boundary: -2; valid label starts from 0!
    vector<pair<double, int> > sorted_vertices(totVertices, pair<double, int>(0, 0));
    // By default the sort function sorts the vector elements on basis of first element of pairs.

    double maxD, minD; // density value Range
    maxD = minD = 0;
    for (int i = 0; i < totVertices; i++) {
        sorted_vertices[i] = pair<double, int>(density[i], i);
        maxD = maxD < density[i] ? density[i] : maxD;
        minD = minD > density[i] ? density[i] : minD;
    }
    std::sort(sorted_vertices.begin(), sorted_vertices.end());
    // By default the sort function sorts the vector elements on basis of first element of pairs.
    // vertices are sorted by the density value ( in fact, it is the inverse height value)

    cout << "Density values: " << maxD << " " << minD << endl;

    // label points by simulated immersion based watershed
    int new_labels = 0; // valid label starts from 0. it must be from 0
    for (auto v : sorted_vertices) {
        int i = v.second; // current vid
        int count = 0;
        int steepest = -1; // vid with the lowest value
        std::set<int> adj_lbls;
        for (auto nv : VV[i]) {
            if (labels[nv] > -1) {
                // adjacent vertex is already labeled
                count++;
                steepest = nv;
                // Note: the steepest is set as the last labeled neighbor, it may not be the steepest
                // however, steepest is only used when there is only one label.
                // the count is the # of labeled neighbors, and it is should be the # of different labels.
                // so use set<int> to collect adjacent labels
                adj_lbls.insert(labels[nv]);
            }
        }
        // label
        count = adj_lbls.size();
        if (count == 0) {
            // no adjacent vertex is labeled -> root of the new region
            regionToMinima[new_labels] = i;
            labels[i] = new_labels;
            new_labels = new_labels + 1;
        } else if (count == 1) {
            // labeled adjacent with the same label -> assign the current vertex with the same label
            labels[i] = labels[steepest];
        } else {
            // labeled adjacent vertices with multiple labels -> the current vertex is boundary
            labels[i] = -2;
        }
    } // finish labeling points

    cout << "Total number of regions: " << new_labels << endl;
    // save the labeled points before the merge
    if (xx_debug) {
        // output_labeled_pts(points, labels, "tree_cluster_points.pts");
    }

    /*
     * Ciccio's old code:
     * persistence value: vid1, vid2
     * persistence from low to high
     * priority_queue<pair<double, pair<int, int>>, vector<pair<double, pair<int, int>>>,std::greater<pair<double,
     * pair<int, int>>>> queue;
     */
    // persistence values between 2 vertices. default is very large value.
    vector<vector<double> > matrixPers(new_labels, vector<double>(new_labels, 10000000));
    // xx: set the queue item to: persistence_value, v1,v2, commonV
    //  add the commonV to correct find the commonV when to merge regions
    priority_queue<pair<double, std::array<int, 3> >, vector<pair<double, std::array<int, 3> > >,
                   std::greater<pair<double, std::array<int, 3> > > >
        queue;

    // persistence based simplification/merge
    for (int i = 0; i < totVertices; i++) {
        if (labels[i] == -2) {
            // boundary point
            double z;
            z = density[i]; // get the value (it is the negative smoothed height) of vertex i
            vector<int> allAdjacent; // adjacent labeled points, (lbl>-1)
            for (auto v : VV[i]) {
                if (labels[v] > -1) {
                    // adjacent vertices is labeled.
                    allAdjacent.push_back(v);
                }
            }
            // check adjacent regions
            for (int v1 = 0; v1 < allAdjacent.size(); v1++) {
                for (int v2 = v1 + 1; v2 < allAdjacent.size(); v2++) {
                    double x1, y1, z1, x2, y2, z2;
                    int l1 = labels[allAdjacent[v1]];
                    int l2 = labels[allAdjacent[v2]];
                    z1 = density[regionToMinima[l1]]; // the value of the root vertex of region
                    z2 = density[regionToMinima[l2]];

                    // FIXME: check the persistence value calculation
                    double xx_per = std::min(std::abs(z - z1), std::abs(z - z2));
                    double ciccio_per = 0.0;
                    if (z1 < z2) {
                        if (matrixPers[l1][l2] > abs(z2 - z)) {
                            matrixPers[l1][l2] = matrixPers[l2][l1] = abs(z - z2);
                            ciccio_per = std::abs(z - z2);
                        }
                    } else {
                        // matrixPers[l1][l2] = matrixPers[l2][l1] > abs(z - z1)
                        if (matrixPers[l1][l2] > abs(z - z1)) {
                            matrixPers[l1][l2] = matrixPers[l2][l1] = abs(z - z1);
                            ciccio_per = std::abs(z - z1);
                        }
                    }
                    // xx_per vs ciccio_per
                    if (xx_debug && (xx_per != ciccio_per)) {
                        cout << z << " " << z1 << " " << z2 << endl;
                        cout << "xx persistence: " << xx_per << " vs ciccio persistence " << ciccio_per << endl;
                    }

                    // ciccio's code
                    // queue.push(pair<double, pair<int, int>>(matrixPers[l1][l2], make_pair(allAdjacent[v1], allAdjacent[v2])));

                    // xx: use xx's persistence calculation
                    matrixPers[l1][l2] = matrixPers[l2][l1] = xx_per;
                    // also add the boundary vertex
                    // https://en.cppreference.com/w/cpp/container/array
                    queue.push(pair<double, std::array<int, 3> >(matrixPers[l1][l2],
                                                                 std::array<int, 3>{
                                                                     allAdjacent[v1], allAdjacent[v2], i
                                                                 }));
                }
            }
        }
    }

    // vids for each valid region (label >-1). lbl: vids.
    // Note: some region may have ZERO vid after the merge process.
    // Note: label must be from 0!
    vector<vector<int> > regionsToVertices(new_labels, vector<int>());
    for (int i = 0; i < totVertices; i++) {
        // boundary points are ignored.
        if (labels[i] > -1) {
            // labeled point
            regionsToVertices[labels[i]].push_back(i);
        }
        if (labels[i] == -1) {
            // not labeled point
            cout << "i: " << i << " is not labeled: " << labels[i] << "\n";
            exit(1);
        }
    }
    // set the persistence value
    double threshold = std::abs(maxD - minD) * perc_pers;
    cout << "maxD, minD, persistence threshold value: " << maxD << " " << minD << " " << threshold << endl;

    while (!queue.empty()) {
        pair<double, std::array<int, 3> > element = queue.top();
        queue.pop();
        int v1 = element.second[0];
        int v2 = element.second[1];
        int xx_commonv = element.second[2];
        int l1 = labels[v1];
        int l2 = labels[v2];
        if (l1 < 0 || l2 < 0) {
            cout << "error! label need to be valid. >=0\n";
            cout << "v1, v2: " << v1 << ", " << v2 << endl;
            cout << "l1, l2: " << l1 << ", " << l2 << endl;
            exit(1);
        }

        if (l1 == l2) // already merge
        {
            continue;
        } else {
            double z, z1, z2;
            z = density[xx_commonv];
            z1 = density[regionToMinima[l1]];
            z2 = density[regionToMinima[l2]];

            // FIXME: check persistence value calculation
            double xx_per = std::min(std::abs(z - z1), std::abs(z - z2));
            double ciccio_per = std::abs(z - std::max(z1, z2));
            if (xx_debug && (xx_per != ciccio_per)) {
                cout << z << " " << z1 << " " << z2 << endl;
                cout << "xx per: " << xx_per << " vs ciccio per" << ciccio_per << endl;
            }
            // xx: use xx's persistence value
            if (element.first != xx_per) {
                // the persistence value need to be updated due to previous merge.
                element.first = xx_per;
                queue.push(element);
            } else {
                if (element.first > threshold) {
                    break;
                    // jump out of checking persistence value
                    // no need to check following ones, as the queue is sorted from low to high persistence values
                }
                // merge regions
                int dominating, dominated; // lower value is the seed (dominating)
                if (z1 > z2) {
                    dominating = l2;
                    dominated = l1;
                } else {
                    dominating = l1;
                    dominated = l2;
                }
                // update labels, and region's vertices
                // cout << labels[v1] << ", " << labels[v2] << " -> ";
                for (auto v : regionsToVertices[dominated]) {
                    labels[v] = dominating;
                }
                // cout << labels[v1] << ", " << labels[v2] << " \n";
                // enlarge dominating region
                regionsToVertices[dominating].insert(regionsToVertices[dominating].begin(),
                                                     regionsToVertices[dominated].begin(),
                                                     regionsToVertices[dominated].end());
                // remove dominated region vertices
                regionsToVertices[dominated].clear();
                // xx: change commonv label as well, remember to add to the regionsToVertices
                // as the commonv is not the boundary anymore
                labels[xx_commonv] = dominating;
                regionsToVertices[dominating].push_back(xx_commonv);
                // cout << "merge: " << dominated << " -> " << dominating << endl;
                new_labels--;
            }
        } // check two adjacent regions
    } // check queue

    cout << "after merge, total number of Tree clusters: " << new_labels << endl;

    // check labeled points -- can be skipped
    if (xx_debug) {
        int tc_num = 0;
        for (int i = 0; i < regionsToVertices.size(); ++i) {
            std::vector<int> t = regionsToVertices[i];
            if (!t.empty()) {
                tc_num = tc_num + 1;
            }
            for (const auto& vid : t) {
                int l = labels[vid];
                if (l != i) {
                    cout << "label is not matched\n";
                    cout << "vid: " << vid << ", lbl: " << l << ", in the region lbl: " << i << endl;
                    exit(1);
                }
            }
        }
        std::set<int> unique_valid_lbls;
        for (int i = 0; i < labels.size(); ++i) {
            int l = labels[i];
            if (l > -1) {
                unique_valid_lbls.insert(l);
                std::vector<int> r = regionsToVertices[l];
                if (std::find(r.begin(), r.end(), i) == r.end()) {
                    cout << "i: " << i << ", lbl: " << l << ", not in region\n";
                    exit(1);
                }
            }
        }

        if (tc_num != unique_valid_lbls.size()) {
            cout << "after merge, final survived labels\n";
            for (const auto& l : unique_valid_lbls) {
                cout << l << " ";
            }
            cout << endl;
            cout << "# of clusters: " << tc_num << ", vs " << unique_valid_lbls.size() << endl;
            cout << "Errors in the labeled points, need to check!\n";
            exit(1);
        }
    } // check merged results. can be skipped

    // update the boundary points.
    // xx: keep the boundary points labeled as boundary points.
    // boundary points will be discarded (not used) in the further cluster segmentation
    // we only change boundary label if all adjacent points with the same label
    for (int i = 0; i < totVertices; i++) {
        if (labels[i] == -2) {
            std::set<int> adj_lbls;
            int tmp_lbl = -1;
            for (auto v : VV[i]) {
                if (labels[v] > -1) {
                    adj_lbls.insert(labels[v]);
                    tmp_lbl = labels[v];
                }
            }
            // xx: only when there is one adjacent label, the boundary point label is changed
            if (adj_lbls.size() == 1) {
                labels[i] = tmp_lbl;
                // xx: no need to update regionsToVertices, as regionsToVertices is no longer needed
            }
        }
    }

    // xx's update: re-organize label
    pts_lbls = vector<int>(labels.size(), -1); // init segment, which is the return value
    std::map<int, int> old2new; // old label to new label
    new_labels = 0; // valid label starts from 0
    for (int i = 0; i < labels.size(); ++i) {
        int old_l = labels[i];
        if (old_l == -1) {
            cout << "error! point: " << i << " is not labeled\n";
            exit(1);
        }
        if (old_l == -2) {
            // boundary. keep it
            pts_lbls[i] = -2;
        } else {
            if (old2new.find(old_l) == old2new.end()) {
                old2new[old_l] = new_labels;
                new_labels = new_labels + 1;
            }
            pts_lbls[i] = old2new[old_l];
        }
    }

    cout << "# of final tree clusters: " << new_labels << endl;

    // have tree clusters: tree cluster id: vids
    for (int i = 0; i < pts_lbls.size(); i++) {
        if (pts_lbls[i] < 0) {
            // valid label starts from 0
            continue;
        }
        if (tree_clusters.find(pts_lbls[i]) == tree_clusters.end()) {
            tree_clusters[pts_lbls[i]] = vector<int>();
        }
        tree_clusters[pts_lbls[i]].push_back(i);
    }

    // note: boundary points will be ignored in the following segmentation
    return 0;
}

int compute_cylinder(double& th_height, const double& radius, const vector<vector<double> >& pts,
                     const vector<int>& pids, vector<vector<int> >& VV, double& zmin, double& zmax) {
    /*
     * connect 3D points to generate VV relationship
     * for each point, the search cylinder is build with the radius. extrude the cylinder from the point. (direction up)
     * the point is the original point at the bottom circle
     * the cylinder height is set from height.
     *
     * Inputs:
     *  pts: all points
     *  pids: point ids in the cluster. global point indices.
     *  height < 0 auto set, otherwise, use height. it is height rate
     * Outputs:
     * VV: index of point indexes in the pids. local vids. thus, a true pid = pids[VV...] todo: check it?
     */

    // here we start working for each cluster
    double minz, maxz, miny, maxy, minx, maxx;
    minz = maxz = miny = maxy = minx = maxx = 0;

    int totVertices = pids.size(); // vertices.size();
    VV = vector<vector<int> >(totVertices, vector<int>());
    map<double, vector<int> > sort_by_height; // height: vids
    // Note: we can ensure the height is saved from low to high?
    // https://www.cplusplus.com/reference/map/map/
    // Internally, the elements in a map are always sorted by
    // its key following a specific strict weak ordering criterion
    // indicated by its internal comparison object (of type Compare).
    for (int i = 0; i < totVertices; i++) {
        double x, y, z;
        x = pts[pids[i]][0];
        y = pts[pids[i]][1];
        z = pts[pids[i]][2];
        if (sort_by_height.find(z) == sort_by_height.end()) {
            // new height found
            sort_by_height[z] = vector<int>();
        }
        sort_by_height[z].push_back(i);
        minz = minz > z ? z : minz;
        maxz = maxz < z ? z : maxz;
        miny = miny > y ? y : miny;
        maxy = maxy < y ? y : maxy;
        minx = minx > x ? x : minx;
        maxx = maxx < x ? x : maxx;
    }
    //    cout << "Interval of values x " << minx << " " << maxx << endl;
    //    cout << "Interval of values y " << miny << " " << maxy << endl;
    cout << "Interval of values z " << minz << " " << maxz << endl;
    zmax = maxz;
    zmin = minz;
    vector<double> heights = vector<double>(sort_by_height.size()); // unique heights
    vector<vector<int> > vertices_h =
        vector<vector<int> >(sort_by_height.size()); // the vertex ids. in the same order with heights
    int i = 0;
    for (auto val : sort_by_height) {
        // height: local vids
        heights[i] = val.first;
        vertices_h[i++] = val.second;
    }

    // compute the height threshold value
    double height = (maxz - minz) * th_height;
    //  if (height <= 0) {
    //    height = (maxz - minz) / 2.0;
    //  }
    cout << "height rate: " << th_height << ". cylinder height: " << height << endl;

    // connect VV
    for (int index = 0; index < heights.size(); index++) {
        for (auto v : vertices_h[index]) {
            // note: here are local vids
            for (int up = index; up < heights.size(); up++) {
                for (auto v2 : vertices_h[up]) {
                    // note: here are local vids
                    if (v != v2) {
                        double x, y, z, x2, y2, z2;
                        x = pts[pids[v]][0];
                        y = pts[pids[v]][1];
                        z = pts[pids[v]][2];
                        x2 = pts[pids[v2]][0];
                        y2 = pts[pids[v2]][1];
                        z2 = pts[pids[v2]][2];
                        if (z < z2) {
                            // if (sqrt(pow(x - x2, 2) + pow(y - y2, 2)) < radius) {
                            if (pow(x - x2, 2) + pow(y - y2, 2) < pow(radius, 2)) {
                                // here check if the points are close enough in xy-plane
                                if (std::abs(z - z2) < std::abs(height)) {
                                    // here check if the points are close enough
                                    // if (pow(z - z2, 2) < pow(height, 2)) {  // here check if the points are close enough
                                    VV[v].push_back(v2);
                                    VV[v2].push_back(v);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int segment_cluster(const vector<vector<double> >& pts, const double& zmin, const double& zmax,
                    const vector<int>& vertices, const vector<vector<int> >& VV, const bool& height_funtype,
                    const double& perc_pers, vector<int>& pts_lbls, int& regionN) {
    /*
     * given the graph of 3D point of cluster points here, segment results
     * Inputs:
     * pts: entire points.
     * vertices: tree cluster points (vids are global indices)
     * VV: graph of tree clusters. vid is the local indices
     * height_funtype:
     *   true: height function. tree bottom focused;
     *   false: reverse height function. treetops focused
     *
     * pts_lbls: entire points label. global vids.
     * perc_pers. real persistence value =  abs(maxZ - minZ) * perc_pers;
     * regionN: the starting region number/id
     *
     * note:
     * VV: uses local vids, based on vertices.
     * update results:
     * segmentation. -> labels of all points (not only the cluster points here)
     */

    // todo: check the local pid and global pids
    int totVertices = vertices.size();
    vector<int> labels(totVertices, -1); // labels of cluster points. default: not labeled -1. local pids.
    map<int, int> regionToMinima; // regioin id: minima pid. local pids.
    vector<pair<double, int> > sorted_vertices(totVertices, pair<double, int>(0, 0));

    for (int i = 0; i < totVertices; i++) {
        double x, y, z;
        x = pts[vertices[i]][0];
        y = pts[vertices[i]][1];
        z = pts[vertices[i]][2];
        if (height_funtype) {
            // height function
            sorted_vertices[i] = pair<double, int>(z, i);
        } else {
            // inverse height function
            sorted_vertices[i] = pair<double, int>(-z, i);
        }
    }
    std::sort(sorted_vertices.begin(), sorted_vertices.end());

    int new_labels = 0; // label starts from 0.
    for (auto v : sorted_vertices) {
        int i = v.second;
        int count = 0;
        int steepest = -1;
        std::set<int> adj_lbls;
        for (auto nv : VV[i]) {
            if (labels[nv] > -1) {
                count++;
                steepest = nv;
                adj_lbls.insert(labels[nv]);
            }
        }
        count = adj_lbls.size();
        if (count == 0) {
            // new minima
            regionToMinima[new_labels] = i;
            labels[i] = new_labels;
            new_labels = new_labels + 1;
        } else if (count == 1) {
            // adjacent to one label
            labels[i] = labels[steepest];
        } else {
            // boundary
            labels[i] = -2;
        }
    }
    if (xx_debug) {
        cout << "Total number of trees " << new_labels << endl;
    }

    vector<vector<double> > matrixPers(new_labels, vector<double>(new_labels, 10000));
    /* ciccio's code:
     * priority_queue<pair<double, pair<int, int>>, vector<pair<double, pair<int, int>>>,
     * std::greater<pair<double, pair<int, int>>>> queue;
     */
    // xx's fix: persistence value: {v1, v2, commonv}
    priority_queue<pair<double, std::array<int, 3> >, vector<pair<double, std::array<int, 3> > >,
                   std::greater<pair<double, std::array<int, 3> > > >
        queue;

    for (int i = 0; i < totVertices; i++) {
        if (labels[i] == -2) {
            // boundary
            vector<int> allAdjacent;
            double x, y, z;
            x = pts[vertices[i]][0];
            y = pts[vertices[i]][1];
            z = pts[vertices[i]][2];
            if (!height_funtype) {
                z = -z;
            }
            for (auto v : VV[i]) {
                // note: vv has local pids
                if (labels[v] > -1) {
                    // points are labeled
                    allAdjacent.push_back(v);
                }
            }
            for (int v1 = 0; v1 < allAdjacent.size(); v1++) {
                for (int v2 = v1 + 1; v2 < allAdjacent.size(); v2++) {
                    double x1, y1, z1, x2, y2, z2;
                    int l1 = labels[allAdjacent[v1]];
                    int l2 = labels[allAdjacent[v2]];
                    x1 = pts[vertices[regionToMinima[l1]]][0];
                    y1 = pts[vertices[regionToMinima[l1]]][1];
                    z1 = pts[vertices[regionToMinima[l1]]][2];
                    x2 = pts[vertices[regionToMinima[l2]]][0];
                    y2 = pts[vertices[regionToMinima[l2]]][1];
                    z2 = pts[vertices[regionToMinima[l2]]][2];

                    if (!height_funtype) {
                        z1 = -z1;
                        z2 = -z2;
                    }

                    /* ciccio's code
                     * if (z1 < z2) {
                     * if (matrixPers[l1][l2] > abs(z2 - z)) matrixPers[l1][l2] = matrixPers[l2][l1] = abs(z2 - z);
                     * } else {
                     * if (matrixPers[l1][l2] > abs(z1 - z)) matrixPers[l1][l2] = matrixPers[l2][l1] = abs(z1 - z);
                     * }
                     */
                    // xx: set persistence value
                    double xx_per = std::min(std::abs(z - z1), std::abs(z - z2));
                    matrixPers[l1][l2] = matrixPers[l2][l1] = xx_per;
                    queue.push(pair<double, std::array<int, 3> >(matrixPers[l1][l2],
                                                                 std::array<int, 3>{
                                                                     allAdjacent[v1], allAdjacent[v2], i
                                                                 }));
                }
            }
        }
    }

    vector<vector<int> > regionsToVertices(new_labels, vector<int>());
    for (int i = 0; i < totVertices; i++) {
        if (labels[i] > -1) {
            regionsToVertices[labels[i]].push_back(i);
        }
    }
    // merge
    double threshold = std::abs(zmax - zmin) * perc_pers;
    cout << "minz, maxz, th_pers: " << zmin << " " << zmax << " " << threshold << endl;
    while (!queue.empty()) {
        pair<double, std::array<int, 3> > element = queue.top();
        queue.pop();
        int v1 = element.second[0];
        int v2 = element.second[1];
        int xx_commonv = element.second[2];
        int l1 = labels[v1];
        int l2 = labels[v2];
        if (l1 == l2) {
            continue;
        } else {
            int commonv;
            double x, y, z, x1, y1, z1, x2, y2, z2;
            x = pts[vertices[xx_commonv]][0];
            y = pts[vertices[xx_commonv]][1];
            z = pts[vertices[xx_commonv]][2];
            x1 = pts[vertices[regionToMinima[l1]]][0];
            y1 = pts[vertices[regionToMinima[l1]]][1];
            z1 = pts[vertices[regionToMinima[l1]]][2];
            x2 = pts[vertices[regionToMinima[l2]]][0];
            y2 = pts[vertices[regionToMinima[l2]]][1];
            z2 = pts[vertices[regionToMinima[l2]]][2];
            if (!height_funtype) {
                z = -z;
                z1 = -z1;
                z2 = -z2;
            }

            double xx_per = std::min(std::abs(z - z1), std::abs(z - z2));
            if (element.first != xx_per) {
                // abs(z - max(z1, z2))
                element.first = xx_per; // abs(z - max(z1, z2));
                queue.push(element);
            } else {
                if (element.first > threshold) {
                    continue;
                }
                int dominating, dominated;
                if (z1 > z2) {
                    dominating = l2;
                    dominated = l1;
                } else {
                    dominating = l1;
                    dominated = l2;
                }
                for (auto v : regionsToVertices[dominated]) {
                    labels[v] = dominating;
                }
                regionsToVertices[dominating].insert(regionsToVertices[dominating].begin(),
                                                     regionsToVertices[dominated].begin(),
                                                     regionsToVertices[dominated].end());
                regionsToVertices[dominated].clear();
                // xx: change commonv label as well, remember to add to the regions2Vertices
                // as the commonv is not the boundary any more
                labels[xx_commonv] = dominating;
                regionsToVertices[dominating].push_back(xx_commonv);

                new_labels--;
            }
        } // check different adjacent regions
    } // check queue

    if (xx_debug) {
        cout << "after merge, total number of trees: " << new_labels << endl;
    }

    // update boundary points
    // ciccio's code
    //  for (int i = 0; i < totVertices; i++) {
    //    if (labels[i] == -2) {
    //      for (auto v : VV[i])
    //        if (labels[v] > -1) {
    //          labels[i] = labels[v];
    //        }
    //    }
    //  }
    // xx
    for (int i = 0; i < totVertices; i++) {
        if (labels[i] == -2) {
            std::set<int> adj_lbls;
            int tmp_lbl = -1;
            for (auto v : VV[i]) {
                if (labels[v] > -1) {
                    adj_lbls.insert(labels[v]);
                    tmp_lbl = labels[v];
                }
            }
            // xx's fix: only when there is one adjacent label, the boundary point label is changed
            if (adj_lbls.size() == 1) {
                labels[i] = tmp_lbl;
            }
        }
    }

    // xx's update: keep boundary points (as label these points -2)
    int newIRegion = regionN;
    map<int, int> newLabels; // old label to new label
    for (auto l : labels) {
        if (l == -2) {
            // keep boundary
            newLabels[l] = -2;
        } else if (newLabels.find(l) == newLabels.end()) {
            newLabels[l] = newIRegion;
            newIRegion = newIRegion + 1;
        }
    }
    regionN = newIRegion;

    for (int i = 0; i < totVertices; i++) {
        pts_lbls[vertices[i]] = newLabels[labels[i]];
    }
    return 0;
}

int merge_labels(const vector<int>& lbls1, const vector<int>& lbls2, vector<int>& lbls, int& trees_num) {
    /*
     * lbls: point labels.
     */
    // here we compute the intersecting segmentation
    // cout<<"we start intersecting\n";
    map<pair<int, int>, int> ms_segm; // <seg label1, seg label2>: final seg id
    int count = 0;
    int pts_num = std::max(lbls1.size(), lbls2.size());
    lbls = vector<int>(pts_num, -1); //  point id:  tree id

    cout << "height + inverse height function\n";
    for (int i = 0; i < lbls1.size(); i++) {
        if (lbls1[i] < 0 || lbls2[i] < 0) {
            lbls[i] = -1; // any boundary points are not considered
            continue;
        }
        pair<int, int> segm = make_pair(lbls1[i], lbls2[i]);
        if (ms_segm.find(segm) == ms_segm.end()) {
            ms_segm[segm] = count++;
        }
        lbls[i] = ms_segm[segm];
    }

    trees_num = count;
    return 1;
}

int post_process(vector<int>& pts_lbls, const int& trees_num, const string& ds = "newfor") {
    /* remove tiny segmentation results
     *
     */

    // sort final labels from # large to small
    vector<int> total_regions(trees_num, 0); // region id: # of pts
    int pts_num = pts_lbls.size();
    for (int i = 0; i < pts_num; i++) {
        if (pts_lbls[i] < 0) {
            continue;
        }
        total_regions[pts_lbls[i]] += 1;
    }

    double avgSize = 0;
    double maxSize = 0;
    std::vector<int> treeSizes;
    for (int i = 0; i < total_regions.size(); i++) {
        treeSizes.emplace_back(total_regions[i]);
        avgSize += total_regions[i];
        if (total_regions[i] > maxSize) {
            maxSize = total_regions[i];
        }
    }
    avgSize /= (float)total_regions.size();
    if (xx_debug) {
        cout << "tree pts average size: " << avgSize << " and max size: " << maxSize << endl;
    }

    // small tree points clouds are marked as -1
    // sort using a lambda expression
    std::sort(treeSizes.begin(), treeSizes.end(), [](int a, int b) { return a > b; });

    int minSize = treeSizes.back();
    // NOTE: may need to adjust thSize for different datasets.
    // NEWFOR: avgSize. default
    // SERC: avgSize and avgSize
    double thSize = 0;

    int xx_count = 0;
    for (auto ts : treeSizes) {
        if (ts >= avgSize) {
            thSize = thSize + ts;
            xx_count++;
        } else {
            break;
        }
    }

    if (ds == "newfor") {
        // for NEWFOR
        thSize = avgSize;
    } else {
        // for SERC
        thSize = thSize / xx_count;
        thSize = thSize * 0.85;
        // thSize = minSize + (maxSize - minSize) * 0.75;
        // thSize = avgSize - (maxSize - avgSize) * 0.01;
    }

    cout << "min size: " << minSize << ",  max size: " << maxSize << endl;
    cout << "average size: " << avgSize << "\n";
    cout << "thSize: " << thSize << endl;

    for (int i = 0; i < pts_num; i++) {
        if (pts_lbls[i] < 0) {
            continue;
        }
        if ((double)total_regions[pts_lbls[i]] < thSize) {
            // avgSize
            pts_lbls[i] = -1;
        }
    }

    vector<int> real_final(pts_num, -1); // new labeled id
    map<int, int> newindex; // old label id : new label id
    map<int, int> count_all; // label: # pts
    int count = 0;
    for (int i = 0; i < pts_num; i++) {
        if (pts_lbls[i] != -1) {
            if (newindex.find(pts_lbls[i]) == newindex.end()) {
                newindex[pts_lbls[i]] = count;
                count_all[count] = 0;
                count++;
            }
            real_final[i] = newindex[pts_lbls[i]];
            count_all[real_final[i]] += 1;
        } else {
            real_final[i] = -1;
        }
    }

    cout << "final trees#: " << count << endl;

    pts_lbls = real_final;
    return 1;
}

void output_labeled_pts(vector<vector<double> >& points, vector<int>& pts_lbls, const string& outfile) {
    /*
     * output labeled points
     * x y z lbl1 lbl2
     *
     * label:
     * -1: cluster boundary
     * -2: boundary of segments inside cluster
     * >=0: valid labels
     */

    cout << "output: " << outfile << endl;

    ofstream ofs(outfile);
    ofs << std::fixed << std::setprecision(2);
    int pts_num = points.size();
    for (int i = 0; i < pts_num; ++i) {
        double x, y, z;
        x = points[i][0];
        y = points[i][1];
        z = points[i][2];
        int l = pts_lbls[i];
        ofs << x << " " << y << " " << z << " " << l << endl;
    }
    ofs.close();
}