#include "formangradient.h"
#include <filesystem> // C++ 17
#include <omp.h>


FormanGradient::FormanGradient(const string &infile, const int &funID) {
    // get infile infos
    _infile = infile;
    std::filesystem::path absolutePath = std::filesystem::absolute(infile);
    _file_name = absolutePath.stem().string();
    std::filesystem::path parentPath = absolutePath.parent_path();
    _workspace = parentPath.string() + "/";
    _file_extension = absolutePath.extension();

    cout << "Infile infos:\n";
    cout << "path: " << absolutePath << endl;
    cout << "file name: " << _file_name << endl;
    cout << "extension: " << _file_extension << endl;

    cout << "\nReading input\n";
    IO_Timer time;
    time.start();
    sc = SimplicialComplex();

    cout << "Entrato" << endl;

    string fileMesh = infile;

    if (fileMesh.find(".off") != string::npos) {
        sc.readOFF(fileMesh.c_str(), funID);
    } else if (
        fileMesh.find(".ply") != string::npos) {
        sc.readPLY(fileMesh.c_str());
    } else {
        cout << "Unknown format" << endl;
        exit(1);
    }
    time.stop();
    cout << "Reading time " << time.getElapsedTime() << " s" << endl;

    time.start();

    cout << "Computing fields" << endl;
    // buildFiltrations(sc);
    // buildFiltrations2(argc,argv);

    int nField = sc.getVertex(0).getCoordinates().size() - 3; // remove xyz
    cout << "Simplicial complex dimension " << sc.getComplexDim() << endl;
    cout << "** Number of fields to create scalar function value: " << nField << ". **" << endl;

    // create gradient encoding
    gradient = GradientEncoding(sc);
    // final filtration
    filtration = vector<uint>(sc.getVerticesNum(), 0);
    // simulation of simplicity for each component
    componentBasedFiltration = vector<vector<uint> >(nField, vector<uint>(sc.getVerticesNum(), 0));
    // input scalar values
    scalarValues = vector<vector<float> >(sc.getVerticesNum(), vector<float>(nField, 0));

    // cout << "Start" << endl;

    vector<pair<float, uint> > injectiveF(sc.getVerticesNum());

    for (int i = 0; i < nField; i++) {
        for (int j = 0; j < sc.getVerticesNum(); j++) {
            float val = sc.getVertex(j).getCoordinate(i + 3);
            scalarValues[j][i] = val;
            // cout << val << endl;
            injectiveF[j] = pair<float, uint>(val, j);
        }
        // cout << "Stop" << endl;

        sort(injectiveF.begin(), injectiveF.end(),
             bind(&FormanGradient::filtrComparer, this, boost::placeholders::_1, boost::placeholders::_2));
        int ind = 0;
        for (auto p: injectiveF) {
            componentBasedFiltration[i][p.second] = ind++;
        }
    }

    vector<vector<uint> > buildFiltration = vector<vector<uint> >(sc.getVerticesNum(), vector<uint>(nField, 0));
    for (int i = 0; i < nField; i++) {
        for (int j = 0; j < sc.getVerticesNum(); j++) {
            buildFiltration[j][i] = componentBasedFiltration[i][j];
        }
    }

    for (int i = 0; i < sc.getVerticesNum(); i++) {
        buildFiltration[i].push_back(i);
    }

    sort(buildFiltration.begin(), buildFiltration.end());

    // filtration created
    int ind = 0;
    for (auto vec: buildFiltration) {
        filtration[vec.back()] = ind++;
    }
    time.stop();
    cout << "Filtration time " << time.getElapsedTime() << " s" << endl;
}

FormanGradient::~FormanGradient() {
}

void FormanGradient::computeFormanGradient(bool computeTopsByBatch) {
    cout << "\nCompute Forman Gradient\n";
    IO_Timer time;

    time.start();
    if (computeTopsByBatch) {
        sc.storeFullStar(); // todo: parallelize it?
    }
    time.stop();

    cout << "Tops computed " << time.getElapsedTime() << " s" << endl;

    //map<pair<int, int>, SSet> filtrationAll; // <- seems used to debug
    //auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this, boost::placeholders::_1, boost::placeholders::_2);
    time.start();

#pragma omp parallel for num_threads(6)
    // todo: set num_threads in a better way, from the user?
    for (uint i = 0; i < filtration.size(); i++) {
        //int thread_id = omp_get_thread_num();
        //std::cout << "Iteration " << i << " is executed by thread " << thread_id << std::endl;


        vector<SSet> lwStars; // [{set of simplices}]
        splitVertexLowerStar(i, lwStars);
        // split the lower star of a vertex v according to the sublevelsets of the function
        // to support multi-parameter/fields

        for (auto lw: lwStars) {
            // uncomment here for old software
            homotopy_expansion(lw); // critical cells generated here

            // debug
            // for (auto s: lw) {
            //     // the code more likes a debug. todo: remove
            //     vector<uint> filtr = this->simplexFiltration(s);
            //     pair<int, int> filrNew(filtr[0], filtr[1]);
            //     if (filtrationAll.find(filrNew) == filtrationAll.end()) filtrationAll[filrNew] = SSet(foo);
            //     filtrationAll[filrNew].insert(s);
            // }
        }
    }
    time.stop();
    cout << "Forman gradient computed " << time.getElapsedTime() << " s" << endl;

    // here I was trying to print out everything for rivet
    // fstream fStream("complexRivet.txt", ios::out);
    //    fStream << "bifiltration" << endl;
    //    fStream << "f1" << endl;
    //    fStream << "f2" << endl;
    //    if (fStream.is_open()) {
    //        for (auto fset: filtrationAll) {
    //            for (auto simpl: fset.second) {
    //                vector<int> indexes = *simpl.getVertices();
    //                for (auto v: indexes)
    //                    fStream << v << " ";
    //                fStream << fset.first.first << " " << fset.first.second << endl;
    //            }
    //        }
    //    }
    //    fStream.close();

    // sc.saveVTK("mesh.vtk");

    cout << "The Forman gradient has critical cells: " << endl;
    for (auto lvl: criticalS) {
        cout << "Dim: " << lvl.first << "  #: " << lvl.second.size() << endl;
    }
    cout << endl;
}

void FormanGradient::splitVertexLowerStar(int v, vector<SSet> &lwStars) {
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this, boost::placeholders::_1, boost::placeholders::_2);
    map<vector<uint>, SSet> lws;
    implicitS vert = implicitS(v);

    vector<uint> filtrLvl = simplexFiltration(vert);

    if (lws.find(filtrLvl) == lws.end()) {
        lws[filtrLvl] = SSet(foo);
    }
    lws[filtrLvl].insert(implicitS(v));

    // cout << sc.getComplexDim() << endl;
    for (int i = 1; i <= sc.getComplexDim(); i++) {
        SSet *lw = vertexLowerStar(v, i);
        // cout << lw->size() << endl;

        for (auto s: *lw) {
            // cout << s << endl;
            vector<uint> filtrLvl = simplexFiltration(s);

            //            for(auto v : filtrLvl)
            //                cout << v << " ";
            //            cout << endl;

            if (lws.find(filtrLvl) == lws.end()) {
                lws[filtrLvl] = SSet(foo);
            }
            lws[filtrLvl].insert(s);
        }
        delete lw;
    }

    for (auto l: lws) lwStars.push_back(l.second);
}

void FormanGradient::homotopy_expansion(SSet &simplexes) {
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this, boost::placeholders::_1, boost::placeholders::_2);
    vector<SSet> sdiv(4, SSet(foo));
    uint alive = 0;

    for (auto s: simplexes) {
        sdiv[s.getDim()].insert(s);
        alive++;
    }

    int d = 1;

    while (alive != 0) {
        if (d < sdiv.size()) {
            while (sdiv[d - 1].size() > 0) {
                // cout << "Here again " << d << endl;

                list<implicitS> toRemove;
                for (auto s: sdiv[d]) {
                    // cout << "Next simplex " << s << endl;
                    implicitS nextPair;

                    bool nmPairable = numPairableLowerStar(s, sdiv[d - 1], nextPair) == 1;

                    if (nmPairable) {
                        // cout << "prima di pair con " << nextPair << endl;
                        setPair(s, nextPair);

                        sdiv[d - 1].erase(nextPair);
                        toRemove.push_back(s);
                        alive -= 2;
                    }
                }

                if (!toRemove.empty()) {
                    // cout << "Prima di rimuovere " << toRemove.size() << " " << sdiv[d].size() << endl;
                    //                    for(auto s : sdiv[d]){
                    //                        //cout << s << endl;
                    //                    }

                    for (auto s: toRemove) {
                        // cout << s << endl;
                        sdiv[d].erase(s);
                    }
                    // cout << sdiv[d].size() << endl;
                } else {
                    implicitS critical = *sdiv[d - 1].begin();
                    sdiv[d - 1].erase(critical);
                    alive--;

#pragma omp critical
                    {
                        if (criticalS.find(critical.getDim()) == criticalS.end()) {
                            SSet crit = SSet(foo);
                            criticalS[critical.getDim()] = crit;
                        }
                        criticalS[critical.getDim()].insert(critical);
                    }
                }
            }

            d++;
        } else {
            implicitS critical = *sdiv[d - 1].begin();
            sdiv[d - 1].erase(critical);
            alive--;

#pragma omp critical
            {
                if (criticalS.find(critical.getDim()) == criticalS.end()) {
                    SSet crit = SSet(foo);
                    criticalS[critical.getDim()] = crit;
                }
                criticalS[critical.getDim()].insert(critical);
            }
        }
    }
}

vector<uint> FormanGradient::simplexFiltration(const implicitS &simpl) {
    vector<uint> filtr(scalarValues[0].size(), 0);
    // nfield items, default value 0

    vector<int> vertices = simpl.getConstVertices();
    for (uint i = 0; i < filtr.size(); i++) {
        for (auto v: vertices) {
            if (componentBasedFiltration[i][v] > filtr[i]) filtr[i] = componentBasedFiltration[i][v];
            // componentBasedFiltration[field_id, pt_id]= sorted_index (aka, filtr_val)
            // xx: filtration value of the simple = max "value" of vertices' "values"
        }
    }

    return filtr;
}

vector<float> FormanGradient::simplexScalarValue(const implicitS &simpl) {
    vector<int> vertis = simpl.getConstVertices();
    vector<float> filtr = scalarValues[vertis[0]];

    for (uint i = 0; i < filtr.size(); i++) {
        for (auto v: vertis) {
            if (scalarValues[v][i] > filtr[i]) filtr[i] = scalarValues[v][i];
        }

        //         float val=0;
        //         for(auto v : vertis){
        //             val +=scalarValues[v][i];
        //         }
        //         val = val/(float)filtr.size();
        //         filtr[i]=val;
    }

    return filtr;
}

int FormanGradient::numPairableLowerStar(const implicitS &next, const SSet &sset, implicitS &pair) {
    vector<implicitS> *boundary = sc.boundaryk(next, next.getDim() - 1);

    int num = 0;
    for (auto s: *boundary) {
        if (sset.find(s) != sset.end()) {
            num++;
            pair = s;
        }
    }

    delete boundary;

    return num;
}

bool FormanGradient::isPaired(const implicitS &simpl) {
    vector<int> vertices = simpl.getConstVertices();
    map<uint, vector<explicitS> > tops;

    // retrieve the fan of top incident into the vertices of next
    for (auto v: vertices) {
        vector<explicitS> *topSimpl = sc.topStar(implicitS(v));
        tops[v] = *topSimpl;
        delete topSimpl;
    }

    vector<explicitS> cobBig = tops[vertices[0]];
    for (int i = 1; i < vertices.size(); i++) {
        vector<explicitS> other = tops[vertices[i]];

        vector<explicitS> result(cobBig.size());
        vector<explicitS>::iterator it =
                set_intersection(cobBig.begin(), cobBig.end(), other.begin(), other.end(), result.begin());
        result.resize(it - result.begin());
        cobBig = result;
    }

    if (cobBig.size() == 0) {
        cout << "the coboundary was empty " << endl;
        cout << simpl.getDim() << endl;

        vector<int> verts = simpl.getConstVertices();

        for (auto v: verts) cout << v << " ";
        cout << endl;

        explicitS top = sc.toExplicit(simpl);
        return gradient.isPaired(simpl, top, sc);
    }

    for (auto s: cobBig) {
        if (gradient.isPaired(simpl, s, sc)) return true;
    }

    return false;
}

void FormanGradient::setPair(const implicitS &next, const implicitS &pair) {
    // next has to be bigger than pair
    assert(next.getDim() > pair.getDim());

    vector<explicitS> *tops = sc.topStar(next);

#pragma omp critical
    {
        for (auto s: *tops) {
            gradient.pair(next, pair, s, sc);
        }
    }
}

void FormanGradient::freePair(const implicitS &next, const implicitS &pair) {
    // next is bigger than pair

    vector<explicitS> *tops = sc.topStar(next);
    for (auto s: *tops) {
        gradient.free(next, pair, s, sc);
    }
}

bool FormanGradient::getPair(const implicitS &next, implicitS &pair) {
    vector<explicitS> *tops = sc.topStar(next);

    // original
    for (auto s: *tops) {
        if (gradient.getPair(next, pair, s, sc)) {
            return true;
        }
    }

    return false;
}

// compute the lower star of a vertex already knowing the top simplexes
SSet *FormanGradient::vertexLowerStar(uint vert, uint dimension) {
    vector<explicitS> *top = sc.topStar(explicitS(0, vert));
    //    cout << "Vertex " << vert << endl;
    //    cout << "Size top " << top->size() << " " << dimension <<endl;
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this, boost::placeholders::_1, boost::placeholders::_2);
    SSet *ret = new SSet(foo);

    for (uint i = 0; i < top->size(); i++) {
        if ((*top)[i].getDim() < dimension) continue;

        if ((*top)[i].getDim() == dimension) {
            vector<int> vertices = sc.getTopSimplex((*top)[i]).getVertices();
            sort(vertices.begin(), vertices.end(),
                 bind(&FormanGradient::sortVerticesFiltration, this, boost::placeholders::_1, boost::placeholders::_2));
            // sortVerticesFiltration: vertex with large filtration value first

            if (vertices[0] == vert) {
                //                for(auto v : vertices)
                //                    cout << v << " ";
                //                cout << endl;

                ret->insert(sc.toImplicit((*top)[i]));
            }
            continue;
        }

        vector<implicitS> *bd = sc.boundaryk(sc.toImplicit((*top)[i]), dimension);
        // cout << bd->size() << endl;

        for (vector<implicitS>::iterator it = bd->begin(); it != bd->end(); it++) {
            vector<int> vertices = it->getConstVertices();
            sort(vertices.begin(), vertices.end(),
                 bind(&FormanGradient::sortVerticesFiltration, this, boost::placeholders::_1, boost::placeholders::_2));

            if (vertices[0] == vert) {
                //                for(auto v : vertices)
                //                    cout << v << " ";
                //                cout << endl;

                ret->insert(*it);
            }
        }

        delete bd;
    }

    delete top;

    return ret;
}


void FormanGradient::saddle2minima(implicitS const &saddle, vector<implicitS> &minima) {
    assert(saddle.getConstVertices().size() == 2); // ensure cell is an edge

    for (int i = 0; i < 2; ++i) {
        // cell.getConstVertices().size()
        int cv = saddle.getConstVertices()[i];
        implicitS s(cv); // vertex simplex
        implicitS next;
        // start from its vertex
        if (getPair(s, next)) {
            assert(isPaired(next) && isPaired(s)); // double check

            stack<implicitS> st_pairs;
            st_pairs.push(next);
            int cur_v = cv;
            while (!st_pairs.empty()) {
                implicitS top = st_pairs.top();
                st_pairs.pop();
                assert(top.getConstVertices().size() == 2); // ensure it is an edge (saddle)
                cur_v = top.getConstVertices()[0] == cur_v ? top.getConstVertices()[1] : top.getConstVertices()[0];
                implicitS ns(cur_v);
                implicitS nnext;
                if (getPair(ns, nnext)) {
                    assert(isPaired(nnext) && isPaired(ns)); // double check
                    assert(nnext.getDim() == 1 && !sc.theSame(nnext, top)); // 1 is from stCell.getDim()
                    st_pairs.push(nnext);
                } else {
                    minima.push_back(ns);
                }
            }
        } else {
            // critical 0-cell is on the 1-cell edge
            minima.push_back(s);
        }
    }

    // note: the above process is so complicated because the implicitS class used
    // default assignment operator, which is memberwise pair aka shallow copy
    // otherwise, iterate way could be better
    // todo: use the iterate, considering with points?
}

void FormanGradient::saddle2maxima(const implicitS &saddle, vector<implicitS> &maxima) {
    // 1 saddle -> 2 max.
    // starting from saddles, --> maxs
    assert(saddle.getConstVertices().size() == 2);
    stack<implicitS> qu; // should be the saddle
    qu.push(saddle);
    //std::vector<implicitS> vpath;
    while (!qu.empty()) {
        implicitS top = qu.top();
        qu.pop();
        vector<implicitS> *bd = sc.coboundaryk(top, top.getDim() + 1); // get triangles. should be 2 triangles, at most
        for (const auto &s: *bd) {
            //cout << "check bd: " << s << endl;
            implicitS next;
            if (getPair(s, next)) {
                assert(isPaired(next) && isPaired(s));
                if (next.getDim() == saddle.getDim() && !sc.theSame(next, top)) {
                    qu.push(next);
                    //vpath.push_back(s);
                    // cout << "next edge: " << next << endl;
                    // cout << "-> tri: " << s << endl;
                }
            } else {
                // s is critical triangle
                if (criticalS[2].find(s) != criticalS[2].end()) {
                    maxima.push_back(s);
                } else {
                    cout << "find max, not in list: " << s << endl;
                }
                //cout << "x-> " << s << " Critical vs " << (criticalS[2].find(s) != criticalS[2].end()) << endl;
                //                if (!vpath.empty()) {
                //                    vpaths.push_back(vpath);
                //                }
                //                vpath.clear();
            }
        }
        delete bd;
    }
}


void FormanGradient::saddle2min_vpaths(const implicitS &saddle, std::vector<std::vector<int> > &vpaths) {
    for (int i = 0; i < 2; ++i) {
        // cell.getConstVertices().size()
        vector<int> vpath;
        int cv = saddle.getConstVertices()[i];
        vpath.push_back(cv);

        implicitS s(cv); // vertex simplex
        implicitS next;
        // start from its vertex
        if (getPair(s, next)) {
            assert(isPaired(next) && isPaired(s)); // double check?

            stack<implicitS> st_pairs;
            st_pairs.push(next);
            int cur_v = cv;
            while (!st_pairs.empty()) {
                implicitS top = st_pairs.top();
                st_pairs.pop();
                assert(top.getConstVertices().size() == 2); // ensure it is an edge (saddle)
                cur_v = top.getConstVertices()[0] == cur_v ? top.getConstVertices()[1] : top.getConstVertices()[0];
                vpath.push_back(cur_v);
                implicitS ns(cur_v);
                implicitS nnext;
                if (getPair(ns, nnext)) {
                    assert(isPaired(nnext) && isPaired(ns)); // double check
                    assert(nnext.getDim() == 1 &&
                        !sc.theSame(nnext, top)); // stCell.getDim() = 1, means 2 vertices. => edge
                    st_pairs.push(nnext);
                }
                //                    else {
                // we meet the critical minima. no paired edge
                // but, the vertex (cur_v) is already inserted to the vpath
                //                        vpath.push_back(ns.getConstVertices()[0]);
                //                    }
            }
        }
        //            else {
        // critical 0-cell is on the 1-cell edge
        // but, the vertex (cv) is already inserted to the vpath
        //                vpath.push_back(s.getConstVertices()[0]);
        //            }

        vpaths.push_back(vpath);
    }
}

void FormanGradient::saddle2max_vpaths(implicitS const &saddle, std::vector<std::vector<implicitS> > &vpaths) {
    // todo: check this code.
    // vpath: edge-triangle-edge-triangles....

    //cout << "input saddle: " << saddle << endl;

    // starting from saddles, --> maxs
    assert(saddle.getConstVertices().size() == 2);
    stack<implicitS> qu; // should be the saddle
    qu.push(saddle);
    std::vector<implicitS> vpath;
    while (!qu.empty()) {
        implicitS top = qu.top();
        qu.pop();
        vector<implicitS> *bd = sc.coboundaryk(top, top.getDim() + 1); // get triangles. should be 2 triangles, at most
        for (const auto &s: *bd) {
            //cout << "check bd: " << s << endl;
            implicitS next;
            if (getPair(s, next)) {
                assert(isPaired(next) && isPaired(s));
                if (next.getDim() == saddle.getDim() && !sc.theSame(next, top)) {
                    qu.push(next);
                    vpath.push_back(s);
                    // cout << "next edge: " << next << endl;
                    // cout << "-> tri: " << s << endl;
                }
            } else {
                // s is critical triangle
                //cout << "x-> " << s << " Critical vs " << (criticalS[2].find(s) != criticalS[2].end()) << endl;
                if (!vpath.empty()) {
                    vpaths.push_back(vpath);
                }
                vpath.clear();
            }
        }
        delete bd;
    }
}
