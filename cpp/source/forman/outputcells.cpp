#include "formangradient.h"

void FormanGradient::out3cells(const list<SSet>& cells, const string& outfile) {
    print_out(outfile.c_str(), cells, 10, 4);
}

void FormanGradient::out2cells(const list<SSet>& cells, bool desc) {
    if (desc)
        print_out("descending2cells.vtk", cells, 5, 3);
    else
        print_out("ascending2cells.vtk", cells, 5, 3);
}

void FormanGradient::out1cells(const list<SSet>& cells, bool desc) {
    if (desc)
        print_out("descending1cells.vtk", cells, 3, 2);
    else {
        print_out("ascending1cells.vtk", cells, 3, 2);
        accurate_asc1cells(cells);
    }
}

void FormanGradient::out0cells(const list<SSet>& cells, const string& filename) {
    print_out((filename + "ascending0cells.vtk").c_str(), cells, 1, 1);
}

void FormanGradient::outCriticalPoints(const SSet& cells) {
    FILE* plotlyFile = fopen("critVal.txt", "w");
    list<Vertex> coords;
    for (auto s : cells) {
        vector<int> vertices = s.getConstVertices();
        Vertex v = sc.getVertex(vertices[0]);

        vector<float> filtr = simplexScalarValue(s);
        fprintf(plotlyFile, "%f %f\n", filtr[0], filtr[1]);

        for (int i = 0; i < vertices.size(); i++)
            v.middlePoint(sc.getVertex(vertices[i]));
        coords.push_back(v);
    }
    fclose(plotlyFile);

    FILE* file = fopen("criticalpoints.vtk", "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", coords.size());

    for (auto c : coords) {
        vector<float> coords = c.getCoordinates();
        // re-scale points in the vtk file
        fprintf(file, "%f %f %f \n", coords[0] - sc.sc_min_x, coords[1] - sc.sc_min_y, coords[2]);
    }

    fprintf(file, "\n\n");

    fprintf(file, "POINT_DATA %d \n", coords.size());
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "cp_type 1 %d float\n", coords.size());

    for (auto s : cells) {
        fprintf(file, "%d ", s.getDim());
    }
    fprintf(file, "\n\n");

    fclose(file);
}

/// Output critical points limited with 0-cells and 1-cells.
/// \param cells: critical simplexes/cells
void FormanGradient::outCriticalPoints_01(const SSet& cells) {
    list<Vertex> coords;
    list<int> cp_type;
    for (const auto& s : cells) {
        if (s.getDim() < 2) {
            vector<int> vertices = s.getConstVertices();
            Vertex v = sc.getVertex(vertices[0]);

            for (int vertice : vertices) {
                v.middlePoint(sc.getVertex(vertice));
            }
            coords.push_back(v);
            cp_type.push_back(s.getDim());
        }
    }

    FILE* file = fopen("criticalpoints_01.vtk", "w");
    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", coords.size());

    for (auto c : coords) {
        vector<float> v_coords = c.getCoordinates();
        fprintf(file, "%f %f %f \n", v_coords[0], v_coords[1], v_coords[2]);
    }

    fprintf(file, "\n\n");

    fprintf(file, "POINT_DATA %d \n", coords.size());
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "cp_type 1 %d float\n", coords.size());

    for (auto s : cp_type) {
        fprintf(file, "%d ", s);
    }
    fprintf(file, "\n\n");

    fclose(file);
}

void FormanGradient::print_out(const char* fileName, list<SSet> const& cells, int param, int dim) {
    map<int, int> oldToNew;
    vector<int> newToOld;

    int triangles = 0;

    int index = 0;
    for (auto sets : cells) {
        triangles += sets.size();
        for (auto s : sets) {
            //cout<<s<<endl;
            vector<int> vertices = s.getConstVertices();
            for (auto v : vertices) {
                if (oldToNew.find(v) == oldToNew.end())
                    oldToNew[v] = index++;
            }
        }
    }

    newToOld = vector<int>(oldToNew.size());
    for (auto v : oldToNew)
        newToOld[v.second] = v.first;

    int vertex_number = oldToNew.size();

    FILE* file = fopen(fileName, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertex_number);

    for (int i = 0; i < vertex_number; i++) {
        vector<float> coords = sc.getVertex(newToOld[i]).getCoordinates();
        fprintf(file, "%f %f %f \n", coords[0], coords[1], coords[2]);
    }

    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", triangles, triangles * (dim + 1));

    for (auto sets : cells) {
        for (auto c : sets) {
            vector<int> vertices = c.getConstVertices();
            fprintf(file, "%d ", dim);
            for (int k = 0; k < dim; k++) {
                fprintf(file, "%d ", oldToNew[vertices[k]]);
            }
            fprintf(file, "\n");
        }
    }
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", triangles);

    for (int i = 0; i < triangles; i++)
        fprintf(file, "%d ", param);
    fprintf(file, "\n\n");

    int nFields = scalarValues[0].size();

    fprintf(file, "POINT_DATA %d \n", vertex_number);
    fprintf(file, "FIELD FieldData %d\n", nFields);

    for (int i = 0; i < nFields; i++) {
        fprintf(file, "originalField_%d 1 %d float\n", i + 1, vertex_number);

        for (int j = 0; j < vertex_number; j++) {
            fprintf(file, "%f ", scalarValues[newToOld[j]][i]);
        }
    }

    fprintf(file, "\n\n");

    fprintf(file, "\n\n");

    fprintf(file, "CELL_DATA %d \n", triangles);
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "descending_2_cells 1 %d int \n", triangles);

    int reg = 0;
    for (auto sets : cells) {
        for (auto c : sets) {
            fprintf(file, "%d ", reg);
        }
        reg++;
    }

    fclose(file);
}

void FormanGradient::accurate_asc1cells(const list<SSet>& cells, const string& vtkfile) {
    cout << "accurate asc1 cells\n";
    map<Vertex, int> vertices;
    list<pair<int, int> > edges; // why list?
    vector<Vertex> vord;

    std::vector<std::vector<pair<int, int> > > gp_es; // grouped edges

    int index = 0;
    for (auto sets : cells) {
        std::vector<pair<int, int> > tmp_es;
        for (auto s : sets) {
            //cout << "\ns " << s << endl;
            vector<explicitS>* tri_top = sc.topStar(s);
            //for (const auto &tt: *tri_top) {
            //    cout << tt << endl;
            //}
            //continue;

            Vertex v1, v2;
            //cout << "size: " << tri_top->size() << endl;
            if (tri_top->size() == 2) {
                v1 = sc.barycenter((*tri_top)[0]);
                v2 = sc.barycenter((*tri_top)[1]);
            } else if (tri_top->size() == 1) {
                v1 = sc.barycenter((*tri_top)[0]);
                v2 = sc.barycenter(s);
            } else {
                cout << "NO " << tri_top->size() << endl;
            }
            delete tri_top;

            //cout << "v1, v2" << v1 << ", " << v2 << endl;
            // continue;

            int indexV1, indexV2;
            // XX: sth wrong. So I make things easy. we don't make sure all center pts are unique.
            // it's not needed for the visualization.
            bool xxD = true;
            if (xxD) {
                vord.push_back(v1);
                vord.push_back(v2);
                edges.push_back(pair<int, int>(vord.size() - 2, vord.size() - 1));

                tmp_es.push_back(pair<int, int>(vord.size() - 2, vord.size() - 1));
                //          vertices[v1] = index;
                //        indexV1 = index;
                //        index=index+1;
                //        vertices[v2] = index;
                //        indexV2 = index;
                //        index=index+1;
                //
                //         cout << "v1 id: " << indexV1 << endl;
                //         cout << "v2 id: " << indexV2 << endl;
                continue;
            } else {
                auto it1 = vertices.find(v1);
                if (it1 == vertices.end()) {
                    vertices[v1] = index;
                    indexV1 = index++;
                } else
                    indexV1 = it1->second;
                cout << "v1 id: " << indexV1 << endl;

                auto it2 = vertices.find(v2);
                if (it2 == vertices.end()) {
                    vertices[v2] = index;
                    indexV2 = index++;
                } else
                    indexV2 = it2->second;
                cout << "v2 id: " << indexV2 << endl;
                edges.push_back(pair<int, int>(indexV1, indexV2));
                cout << "push " << indexV1 << ", " << indexV2 << endl;
            }
        }
        cout << "sets#: " << sets.size() << ", vs tmp es: " << tmp_es.size() << endl;
        gp_es.push_back(tmp_es);
    }

    //exit(1);
    cout << "vertices #: " << vertices.size() << endl;
    cout << "vertices #: " << vord.size() << endl;
    cout << "edges #: " << edges.size() << endl;
    cout << "vpaths #: " << gp_es.size() << endl;
    for (int i = 0; i < 3; ++i) {
        cout << "size: " << gp_es[i].size() << endl;
    }
    //vector<Vertex> vord(vertices.size());
    //for (auto v : vertices) vord[v.second] = v.first;

    string outfile = "ascending1cells_correct.vtk";
    if (!vtkfile.empty()) {
        outfile = vtkfile;
    }
    //FILE *file = fopen("ascending1cells_correct.vtk", "w");
    FILE* file = fopen(outfile.c_str(), "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vord.size());

    for (auto v : vord) {
        fprintf(file, "%f %f %f \n", v.getCoordinate(0), v.getCoordinate(1), v.getCoordinate(2));
    }

    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", edges.size(), edges.size() * 3);

    for (auto e : edges) {
        fprintf(file, "2 %d %d\n", e.first, e.second);
    }
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", edges.size());

    for (int i = 0; i < edges.size(); i++)
        fprintf(file, "3 ");
    fprintf(file, "\n\n");

    fprintf(file, "CELL_DATA %d \n", edges.size());
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "ascending_cells 1 %d int \n", edges.size());

    for (int i = 0; i < gp_es.size(); ++i) {
        for (int j = 0; j < gp_es[i].size(); ++j) {
            fprintf(file, "%d ", i);
        }
    }
    fprintf(file, "\n\n");

    fclose(file);
}

// void FormanGradient::printGraph(MIG& mig){

//    map<Node*,int> nodeToIndex;
//    vector<Vertex> vertices(mig.getNodes().size());

//    cout << "Going to print graph" << endl;

//    int index=0;
//    for(auto n : mig.getNodes()){
//        assert(n!=NULL);
//        CritSet all = n->getCSimplexes();
//        bool first=true;
//        Vertex v;

//        for(auto cs : all){
//            for(auto s : cs){
//                if(first){
//                    v = Vertex(sc.barycenter(s).getCoordinates());
//                    first=false;
//                }
//                else{

//                    Vertex v1 = sc.barycenter(s);
//                    v.middlePoint(v1);
//                }
//            }
//        }

//        nodeToIndex[n]=index;
//        vertices[index++]=v;
//    }

//    cout << "Nodes prepared " << endl;

//    FILE* file = fopen("incidenceGraph.vtk","w");

//    fprintf(file, "# vtk DataFile Version 2.0\n\n");
//    fprintf(file, "ASCII \n");
//    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
//    fprintf(file, "POINTS %d float\n", vertices.size());

//    for(auto v : vertices){
//        fprintf(file, "%f %f %f \n", v.getCoordinate(0),v.getCoordinate(1),v.getCoordinate(2));
//    }

//    fprintf(file, "\n\n");

//    fprintf(file, "CELLS %d %d\n", mig.getArcs().size(), mig.getArcs().size()*3);

//    for(auto e : mig.getArcs()){
//            fprintf(file, "2 %d %d\n", nodeToIndex[e->getNodeUp()], nodeToIndex[e->getNodeDown()]);
//    }
//    fprintf(file, "\n");

//    fprintf(file, "CELL_TYPES %d\n", mig.getArcs().size());

//    for(int i=0; i<mig.getArcs().size(); i++)
//        fprintf(file, "3 ");
//    fprintf(file, "\n\n");

//    fclose(file);

//}

void FormanGradient::saveScalarFieldVTK(const char* filename) {
    // scaled version: points are rescaled.
    FILE* file = fopen(filename, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", sc.getVerticesNum());

    for (int i = 0; i < sc.getVerticesNum(); i++) {
        vector<float> coords = sc.getVertex(i).getCoordinates();
        // re-scale points for the vtk file
        fprintf(file, "%f %f %f \n", coords[0] - sc.sc_min_x, coords[1] - sc.sc_min_y, coords[2]);
    }

    fprintf(file, "\n\n");

    int sizeTop = 0;
    int sizeIndices = 0;
    for (auto d : sc.getTopSimplexesSet()) {
        sizeTop += sc.getTopSimplexesNum(d);
        sizeIndices += sc.getTopSimplexesNum(d) * (d + 2);
    }

    fprintf(file, "CELLS %d %d\n", sizeTop, sizeIndices);

    for (auto d : sc.getTopSimplexesSet()) {
        for (auto s : sc.getTopSimplices(d)) {
            vector<int> verts = s.getVertices();
            fprintf(file, "%d ", verts.size());
            for (int i = 0; i < verts.size(); i++)
                fprintf(file, "%d ", verts[i]);
            fprintf(file, "\n");
        }
    }

    fprintf(file, "CELL_TYPES %d\n", sizeTop);
    for (auto d : sc.getTopSimplexesSet()) {
        for (auto s : sc.getTopSimplices(d)) {
            if (s.getDimension() == 1)
                fprintf(file, "3 ");
            else if (s.getDimension() == 2)
                fprintf(file, "5 ");
            else if (s.getDimension() == 3)
                fprintf(file, "10 ");
        }
    }
    fprintf(file, "\n");

    int nFields = scalarValues[0].size();

    fprintf(file, "POINT_DATA %d \n", sc.getVerticesNum());
    fprintf(file, "FIELD FieldData %d\n", nFields);

    for (int i = 0; i < nFields; i++) {
        fprintf(file, "originalField_%d 1 %d float\n", i + 1, sc.getVerticesNum());

        for (int j = 0; j < sc.getVerticesNum(); j++) {
            fprintf(file, "%f ", scalarValues[j][i]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void FormanGradient::saveScalarFieldVTK_OG(const char* filename) {
    FILE* file = fopen(filename, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", sc.getVerticesNum());

    for (int i = 0; i < sc.getVerticesNum(); i++) {
        vector<float> coords = sc.getVertex(i).getCoordinates();
        // re-scale points for the vtk file
        fprintf(file, "%f %f %f \n", coords[0], coords[1], coords[2]);
    }

    fprintf(file, "\n\n");

    int sizeTop = 0;
    int sizeIndices = 0;
    for (auto d : sc.getTopSimplexesSet()) {
        sizeTop += sc.getTopSimplexesNum(d);
        sizeIndices += sc.getTopSimplexesNum(d) * (d + 2);
    }

    fprintf(file, "CELLS %d %d\n", sizeTop, sizeIndices);

    for (auto d : sc.getTopSimplexesSet()) {
        for (auto s : sc.getTopSimplices(d)) {
            vector<int> verts = s.getVertices();
            fprintf(file, "%d ", verts.size());
            for (int i = 0; i < verts.size(); i++)
                fprintf(file, "%d ", verts[i]);
            fprintf(file, "\n");
        }
    }

    fprintf(file, "CELL_TYPES %d\n", sizeTop);
    for (auto d : sc.getTopSimplexesSet()) {
        for (auto s : sc.getTopSimplices(d)) {
            if (s.getDimension() == 1)
                fprintf(file, "3 ");
            else if (s.getDimension() == 2)
                fprintf(file, "5 ");
            else if (s.getDimension() == 3)
                fprintf(file, "10 ");
        }
    }
    fprintf(file, "\n");

    int nFields = scalarValues[0].size();

    fprintf(file, "POINT_DATA %d \n", sc.getVerticesNum());
    fprintf(file, "FIELD FieldData %d\n", nFields);

    for (int i = 0; i < nFields; i++) {
        fprintf(file, "originalField_%d 1 %d float\n", i + 1, sc.getVerticesNum());

        for (int j = 0; j < sc.getVerticesNum(); j++) {
            fprintf(file, "%f ", scalarValues[j][i]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}