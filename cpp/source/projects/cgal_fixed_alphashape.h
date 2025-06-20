//
// Created by alex on 6/29/21.
//

#ifndef CGAL_AS_CGAL_FIXED_ALPHASHAPE_H
#define CGAL_AS_CGAL_FIXED_ALPHASHAPE_H

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
// fixed alpha shape
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
//
#include <CGAL/Object.h>
#include <limits.h>
#include <stdlib.h>
#include <sys/stat.h>  //mkdir

#include <algorithm>
#include <cassert>
#include <cmath>    /* sqrt */
#include <cstdlib>  // realpath
#include <fstream>
#include <iomanip>  // setprecision
#include <iomanip>  // std::setprecision
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>  // stringstream
#include <sstream>
#include <vector>

// todo: try to make a single header file, remove other .h files.
#include "../usage/xxTimer.h"
#include "../usage/Usage.h"

// todo: reconsider io support in cgal and entire project
// refactor the code
#include "../projects/happly.h"

#include <filesystem>

void recreateFolder(const std::string &folderPath) {
    if (std::filesystem::exists(folderPath)) {
        cout << "found the old folder. delete\n";
        std::filesystem::remove_all(folderPath);
    }
    cout << "create the folder\n";
    std::filesystem::create_directory(folderPath);
}

///
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

// fixed alpha shape
// example:
// https://doc.cgal.org/latest/Alpha_shapes_3/Alpha_shapes_3_2ex_fixed_weighted_alpha_shapes_3_8cpp-example.html#_a4
/// vertex base
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb_Tri;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<K, Vb_Tri> Vb_Fixed;
/// cell base
// https://doc.cgal.org/latest/Alpha_shapes_3/classCGAL_1_1Fixed__alpha__shape__cell__base__3.html
//  default: Delaunay_triangulation_cell_base_3<Traits>
// typedef CGAL::Delaunay_triangulation_cell_base_3<K> Dcb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<K> Cb_Fixed;
/// TDS, and fixed alpha shape
typedef CGAL::Triangulation_data_structure_3<Vb_Fixed, Cb_Fixed> Tds_Fixed;
typedef CGAL::Delaunay_triangulation_3<K, Tds_Fixed, CGAL::Fast_location> Delaunay_Fixed;
typedef CGAL::Fixed_alpha_shape_3<Delaunay_Fixed> Fixed_alpha_shape_3;

using namespace std;

// Fixed Alpha shape

// IO
// read
template<typename T>
void str2vector(const string &str, vector<T> &vec) {
    istringstream iss(str);
    // https://stackoverflow.com/questions/9986091/how-do-you-convert-a-string-into-an-array-of-floats
    std::copy(std::istream_iterator<T>(iss), std::istream_iterator<T>(), std::back_inserter(vec));
}

// output alpha shape in off
void xx_write_fixed_as_off(const Fixed_alpha_shape_3 &as, const std::map<Point, int> &ptlbls,
                           const std::string &as_off_file) {
    // output alpha shape as off file.
    //
https
: // github.com/CGAL/cgal/blob/74af3e7d33e7354e4090a1b17941a3f0c372ca8d/Alpha_shapes_3/examples/Alpha_shapes_3/visible_alpha_shape_facets_to_OFF.cpp
    //  new way to output off
    /*
     * https://doc.cgal.org/latest/Alpha_shapes_3/index.html
     *
     * Types of k-faces of the triangulation: EXTERIOR, SINGULAR, REGULAR or INTERIOR
     * ? what's the boundary?
     *  "touch" the rest world without alpha shape. (just the common sense of boundary)
     * On the boundary                                          | Not on the boundary
     * Regular:                                                 | Interior:
     *   it is a subface of the alpha-complex                   |
     *   which is a subface of a (k+1)-face of the alpha complex|
     * Singular:                                                |
     *   not a subface of a (k+1)-face of the alpha complex     |Exterior: squared radius > given alpha
     */

    // https://doc.cgal.org/latest/Alpha_shapes_3/Alpha_shapes_3_2ex_fixed_weighted_alpha_shapes_3_8cpp-example.html
    std::cout << "Information of the Alpha-Complex:\n";
    std::vector<Fixed_alpha_shape_3::Cell_handle> cells; // tetra
    std::vector<Fixed_alpha_shape_3::Facet> facets; // triangles
    std::vector<Fixed_alpha_shape_3::Edge> edges; // edges

    // tetrahedron = cell, they should be in the interior as it is inside the 3D space
    as.get_alpha_shape_cells(std::back_inserter(cells), Fixed_alpha_shape_3::INTERIOR);

    // triangles
    // tetrahedron will contain regular triangles
    // as regular triangles are included in k+1 simplex (aka tetra) of alpha complex
    // it can be ignored for visualization purpose
    // but, it may be important to output when to analyze the alpha shape
    // as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
    // singular triangles are boundary triangles which are NOT included in any k+1 simplex of the complex
    as.get_alpha_shape_facets(std::back_inserter(facets), Fixed_alpha_shape_3::SINGULAR);

    // edges
    // as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_edges(std::back_inserter(edges), Fixed_alpha_shape_3::SINGULAR);

    std::cout << "The alpha-complex has: " << std::endl;
    std::cout << cells.size() << " interior tetrahedrons" << std::endl;
    std::cout << facets.size() << " singular triangles" << std::endl;
    std::cout << edges.size() << " singular edges" << std::endl;

    // https://github.com/SenthilCaesar/Protein-Complex-Alpha-Shape/blob/master/AlphaShape.cpp
    std::vector<Fixed_alpha_shape_3::Vertex_handle> vertices; //
    as.get_alpha_shape_vertices(std::back_inserter(vertices), Fixed_alpha_shape_3::SINGULAR);
    std::cout << vertices.size() << " singular vertices" << std::endl;

    // vertices: points <-> id
    std::map<Point, size_t> points;
    size_t index = 0;
    // finite_.. is from DT class
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++) {
        points[v_it->point()] = index;
        index++;
    }

    // cout << "# of pts: " << index << " vs " << points.size() << endl;

    size_t tetra_num, tri_num, edge_num;
    tetra_num = cells.size();
    tri_num = facets.size();
    edge_num = edges.size();

    // write file
    std::ofstream ofs(as_off_file);
    ofs << fixed << setprecision(3);
    ofs << "OFF\n" << index << " " << tetra_num + tri_num + edge_num << " 0 \n";
    // write points
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++) {
        Point p = v_it->point();
        int lbl = 0;
        if (ptlbls.find(p) != ptlbls.end()) {
            lbl = ptlbls.at(p);
        }
        //
        if (lbl == 0) {
            ofs << v_it->point() << std::endl;
        } else {
            ofs << p << " " << lbl << std::endl;
        }
    }
    // tetra
    for (const auto &cell: cells) {
        size_t v0 = points.find(cell->vertex(0)->point())->second;
        size_t v1 = points.find(cell->vertex(1)->point())->second;
        size_t v2 = points.find(cell->vertex(2)->point())->second;
        size_t v3 = points.find(cell->vertex(3)->point())->second;
        ofs << "4 " << v0 << " " << v1 << " " << v2 << " " << v3 << std::endl;
    }
    // triangles
    /* https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#ad6a20b45e66dfb690bfcdb8438e9fcae
     * (c,i) is the facet of c opposite to the vertex of index i.
     * c: cell <- tetrahedron; i: vertex
     */
    for (auto &facet: facets) {
        ofs << "3 ";
        auto tmp_tetra = facet.first;
        for (int i = 0; i < 4; i++) {
            if (i != facet.second) {
                ofs << points.find(tmp_tetra->vertex(i)->point())->second << " ";
            }
        }
        ofs << std::endl;
    }
    // edges
    /* https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#af31db7673a6d7d28c0bb90a3115ac695
     * (c,i,j) is the edge of cell c whose vertices indices are i and j.
     */
    for (auto e: edges) {
        ofs << "2 ";
        auto tmp_tetra = e.get<0>();
        int p1, p2;
        p1 = e.get<1>();
        p2 = e.get<2>();
        ofs << points.find(tmp_tetra->vertex(p1)->point())->second << " "
                << points.find(tmp_tetra->vertex(p2)->point())->second << std::endl;
    }

    ofs << std::endl;

    ofs.close();
} // write_as_off


void xx_write_fixed_as_ply_binary(const Fixed_alpha_shape_3 &as,
                                  const std::string &as_out_file) {
    // todo: review the code of extracting simplexes, need to fit the forman gradient calculation
    // output alpha shape as binary ply file
    //
https
: // github.com/CGAL/cgal/blob/74af3e7d33e7354e4090a1b17941a3f0c372ca8d/Alpha_shapes_3/examples/Alpha_shapes_3/visible_alpha_shape_facets_to_OFF.cpp
    //  new way to output off
    /*
     * https://doc.cgal.org/latest/Alpha_shapes_3/index.html
     *
     * Types of k-faces of the triangulation: EXTERIOR, SINGULAR, REGULAR or INTERIOR
     * ? what's the boundary?
     *  "touch" the rest world without alpha shape. (just the common sense of boundary)
     * On the boundary                                          | Not on the boundary
     * Regular:                                                 | Interior:
     *   it is a subface of the alpha-complex                   |
     *   which is a subface of a (k+1)-face of the alpha complex|
     * Singular:                                                |
     *   not a subface of a (k+1)-face of the alpha complex     |Exterior: squared radius > given alpha
     */

    // https://doc.cgal.org/latest/Alpha_shapes_3/Alpha_shapes_3_2ex_fixed_weighted_alpha_shapes_3_8cpp-example.html
    std::cout << "Information of the Alpha-Complex:\n";
    std::vector<Fixed_alpha_shape_3::Cell_handle> cells; // tetra
    std::vector<Fixed_alpha_shape_3::Facet> facets; // triangles
    std::vector<Fixed_alpha_shape_3::Edge> edges; // edges

    // tetrahedron = cell, they should be in the interior as it is inside the 3D space
    as.get_alpha_shape_cells(std::back_inserter(cells), Fixed_alpha_shape_3::INTERIOR);

    // triangles
    // tetrahedron will contain regular triangles
    // as regular triangles are included in k+1 simplex (aka tetra) of alpha complex
    // it can be ignored for visualization purpose
    // but, it may be important to output when to analyze the alpha shape
    // as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
    // singular triangles are boundary triangles which are NOT included in any k+1 simplex of the complex
    as.get_alpha_shape_facets(std::back_inserter(facets), Fixed_alpha_shape_3::SINGULAR);

    // edges
    // as.get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::REGULAR);
    as.get_alpha_shape_edges(std::back_inserter(edges), Fixed_alpha_shape_3::SINGULAR);

    std::cout << "The alpha-complex has: " << std::endl;
    std::cout << cells.size() << " interior tetrahedrons" << std::endl;
    std::cout << facets.size() << " singular triangles" << std::endl;
    std::cout << edges.size() << " singular edges" << std::endl;

    // https://github.com/SenthilCaesar/Protein-Complex-Alpha-Shape/blob/master/AlphaShape.cpp
    std::vector<Fixed_alpha_shape_3::Vertex_handle> vertices; //
    as.get_alpha_shape_vertices(std::back_inserter(vertices), Fixed_alpha_shape_3::SINGULAR);
    std::cout << vertices.size() << " singular vertices" << std::endl;

    // vertices: points <-> id
    std::map<Point, size_t> points;
    size_t index = 0;
    // finite_.. is from DT class
    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++) {
        points[v_it->point()] = index;
        index++;
    }

    // cout << "# of pts: " << index << " vs " << points.size() << endl;
    // size_t tetra_num, tri_num, edge_num;
    // tetra_num = cells.size();
    // tri_num = facets.size();
    // edge_num = edges.size();

    // write file
    // std::vector<std::array<int, 4> > as_tetras;
    // std::vector<std::array<int, 3> > as_triangles;
    // std::vector<std::array<int, 2> > as_edges;
    std::vector<std::array<double, 3> > as_vertices;

std:
    vector<std::vector<int> > ply_total_faces;

    for (auto v_it = as.finite_vertices_begin(); v_it != as.finite_vertices_end(); v_it++) {
        Point p = v_it->point();
        as_vertices.push_back({p.x(), p.y(), p.z()});
    }
    // tetra
    for (const auto &cell: cells) {
        size_t v0 = points.find(cell->vertex(0)->point())->second;
        size_t v1 = points.find(cell->vertex(1)->point())->second;
        size_t v2 = points.find(cell->vertex(2)->point())->second;
        size_t v3 = points.find(cell->vertex(3)->point())->second;
        //as_tetras.push_back({v0, v1, v2, v3});
        ply_total_faces.push_back({v0, v1, v2, v3});
    }
    // triangles
    /* https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#ad6a20b45e66dfb690bfcdb8438e9fcae
     * (c,i) is the facet of c opposite to the vertex of index i.
     * c: cell <- tetrahedron; i: vertex
     */
    for (auto &facet: facets) {
        std::array<int, 3> tmp_face;
        int cnt = 0;
        auto tmp_tetra = facet.first;
        for (int i = 0; i < 4; i++) {
            if (i != facet.second) {
                tmp_face[cnt] = points.find(tmp_tetra->vertex(i)->point())->second;
                cnt = cnt + 1;
            }
        }
        //as_triangles.push_back(tmp_face);
        ply_total_faces.push_back({tmp_face[0], tmp_face[1], tmp_face[2]});
    }
    // edges
    /* https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#af31db7673a6d7d28c0bb90a3115ac695
     * (c,i,j) is the edge of cell c whose vertices indices are i and j.
     */
    for (auto e: edges) {
        auto tmp_tetra = e.get<0>();
        int p1, p2;
        p1 = e.get<1>();
        p2 = e.get<2>();
        size_t v0 = points.find(tmp_tetra->vertex(p1)->point())->second;
        size_t v1 = points.find(tmp_tetra->vertex(p2)->point())->second;
        //as_edges.push_back({v0, v1});
        ply_total_faces.push_back({v0, v1});
    }

    std::cout << "Writing binary PLY: " << as_out_file << std::endl;

    happly::PLYData plyOut;
    plyOut.addVertexPositions(as_vertices);
    plyOut.addFaceIndices(ply_total_faces);


    // Add vertices
    // std::vector<double> x, y, z;
    // for (auto &v: as_vertices) {
    //     x.push_back(v[0]);
    //     y.push_back(v[1]);
    //     z.push_back(v[2]);
    // }
    // plyOut.addElement("vertex", vertices.size());
    // plyOut.getElement("vertex").addProperty<double>("x", x);
    // plyOut.getElement("vertex").addProperty<double>("y", y);
    // plyOut.getElement("vertex").addProperty<double>("z", z);
    //
    // // Add edges
    // if (!as_edges.empty()) {
    //     plyOut.addElement("edge", as_edges.size());
    //     plyOut.getElement("edge").addListProperty<int>("vertex_indices", as_edges);
    // }
    //
    // // Add triangles
    // if (!as_triangles.empty()) {
    //     plyOut.addElement("face", as_triangles.size());
    //     plyOut.getElement("face").addListProperty<int>("vertex_indices", as_triangles);
    // }
    //
    // // Add tetrahedrons
    // if (!as_tetras.empty()) {
    //     plyOut.addElement("tetrahedron", as_tetras.size());
    //     plyOut.getElement("tetrahedron").addListProperty<int>("vertex_indices", as_tetras);
    // }

    // Write to binary file
    plyOut.write(as_out_file, happly::DataFormat::Binary);

    cout << "wrote: " << as_out_file << endl;
} // write_as_ply_binary


// output object from cells, facets and edges in off
void write_object_off(std::vector<Fixed_alpha_shape_3::Cell_handle> cells,
                      std::vector<Fixed_alpha_shape_3::Facet> facets, std::vector<Fixed_alpha_shape_3::Edge> edges,
                      const string &outfile) {
    /*
     * output object based on cells (tetras), facets (triangles), edges (edges)
     */

    // zmin, zmax to remove
    double zmin, zmax;
    zmin = 100000;
    zmax = 0;

    // vertices: points <-> id
    std::map<Point, size_t> pt2id; // point to id
    std::vector<Point> points;
    // tetras
    for (const auto &cell: cells) {
        for (int i = 0; i < 4; ++i) {
            Point p = cell->vertex(i)->point();
            if (pt2id.find(p) == pt2id.end()) {
                points.push_back(p);
                pt2id[p] = points.size() - 1;
                zmin = std::min(zmin, p.z());
                zmax = std::max(zmax, p.z());
            }
        }
    } // end tetra
    // tris
    for (const auto &tri: facets) {
        auto tmp_tetra = tri.first;
        for (int i = 0; i < 4; i++) {
            if (i != tri.second) {
                Point p = tmp_tetra->vertex(i)->point();
                if (pt2id.find(p) == pt2id.end()) {
                    points.push_back(p);
                    pt2id[p] = points.size() - 1;
                    zmin = std::min(zmin, p.z());
                    zmax = std::max(zmax, p.z());
                }
            }
        }
    } // end tris
    // edges
    for (const auto &edge: edges) {
        auto tmp_tetra = edge.get<0>();
        int pid1, pid2;
        pid1 = edge.get<1>();
        pid2 = edge.get<2>();
        Point p1 = tmp_tetra->vertex(pid1)->point();
        Point p2 = tmp_tetra->vertex(pid2)->point();
        if (pt2id.find(p1) == pt2id.end()) {
            points.push_back(p1);
            pt2id[p1] = points.size() - 1;
            zmin = std::min(zmin, p1.z());
            zmax = std::max(zmax, p1.z());
        }
        if (pt2id.find(p2) == pt2id.end()) {
            points.push_back(p2);
            pt2id[p2] = points.size() - 1;
            zmin = std::min(zmin, p2.z());
            zmax = std::max(zmax, p2.z());
        }
    }
    // cout << "z range: " << zmin << " - " << zmax << endl;

    // set output requirements from augments -> output without checking, filtering is done by other code

    //    if (zmin > 3.0) {
    //      cout << "not a valid object\n";
    //      return;
    //    }
    //
    //    if (zmax - zmin < 2) {
    //      cout << "not a valid object\n";
    //      return;
    //    }

    // cout << "# of pts: " << index << " vs " << points.size() << endl;

    cout << "start write: " << outfile << endl;
    size_t tetra_num, tri_num, edge_num;
    tetra_num = cells.size();
    tri_num = facets.size();
    edge_num = edges.size();

    cout << "with:\n";
    cout << "# of tetra: " << tetra_num << endl;
    cout << "# of tris: " << tri_num << endl;
    cout << "# of edges: " << edge_num << endl;

    // write file
    std::ofstream ofs(outfile);
    ofs << fixed << setprecision(3);
    ofs << "OFF\n" << points.size() << " " << tetra_num + tri_num + edge_num << " 0 \n";
    // write points
    for (const auto &p: points) {
        ofs << p << std::endl;
    }
    // tetra
    for (const auto &cell: cells) {
        size_t v0 = pt2id.find(cell->vertex(0)->point())->second;
        size_t v1 = pt2id.find(cell->vertex(1)->point())->second;
        size_t v2 = pt2id.find(cell->vertex(2)->point())->second;
        size_t v3 = pt2id.find(cell->vertex(3)->point())->second;
        ofs << "4 " << v0 << " " << v1 << " " << v2 << " " << v3 << std::endl;
    }
    // triangles
    /* https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#ad6a20b45e66dfb690bfcdb8438e9fcae
     * (c,i) is the facet of c opposite to the vertex of index i.
     * c: cell <- tetrahedron; i: vertex
     */
    for (auto &facet: facets) {
        ofs << "3 ";
        auto tmp_tetra = facet.first;
        for (int i = 0; i < 4; i++) {
            if (i != facet.second) {
                ofs << pt2id.find(tmp_tetra->vertex(i)->point())->second << " ";
            }
        }
        ofs << std::endl;
    }
    // edges
    /* https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#af31db7673a6d7d28c0bb90a3115ac695
     * (c,i,j) is the edge of cell c whose vertices indices are i and j.
     */
    for (auto e: edges) {
        ofs << "2 ";
        auto tmp_tetra = e.get<0>();
        int p1, p2;
        p1 = e.get<1>();
        p2 = e.get<2>();
        ofs << pt2id.find(tmp_tetra->vertex(p1)->point())->second << " "
                << pt2id.find(tmp_tetra->vertex(p2)->point())->second << std::endl;
    }

    ofs << std::endl;

    ofs.close();
}; // write_as_off with cells, facets, and edges

//

//


void get_connected_objects_3(Fixed_alpha_shape_3 &as, std::map<Point, int> &ptlbls, const string &outdir = "",
                             const bool _savelbl = false, const bool _showlog = false) {
    /*
     * idea:
     * instances: cell, facets, edges
     * generate the vertex 2 instances map :)
     *
     * output:
     * 1. labeled points. valid label from 1
     * 2. each connected as files
     *
     */

    float st_memory, end_memory;
    st_memory = xx_getMemoryInstantUsage();
    cout << "Apply get_connected_objects_3\n";

    std::vector<Fixed_alpha_shape_3::Cell_handle> interior_tetras; // tetra
    std::vector<Fixed_alpha_shape_3::Facet> singular_tris; // triangles
    std::vector<Fixed_alpha_shape_3::Edge> singular_edges; // edges

    as.get_alpha_shape_cells(std::back_inserter(interior_tetras), Fixed_alpha_shape_3::INTERIOR);
    as.get_alpha_shape_facets(std::back_inserter(singular_tris), Fixed_alpha_shape_3::SINGULAR);
    as.get_alpha_shape_edges(std::back_inserter(singular_edges), Fixed_alpha_shape_3::SINGULAR);

    int tetras_num = interior_tetras.size();
    int tris_num = singular_tris.size();
    int edges_num = singular_edges.size();

    // label vertex
    size_t count = 0;
    for (Fixed_alpha_shape_3::Finite_vertices_iterator it = as.finite_vertices_begin(); it != as.finite_vertices_end();
         it++) {
        it->info() = count++;
    }

    // vertex to instances (index) based on interior_tetras, singular_tris, singular_edges list
    //  ofstream ofs1("as_instance_pts.txt");
    // ofs1 << fixed << setprecision(3);
    std::map<int, std::set<int> > v2ins; // vertex id to instances ids.
    // std::map<Point, std::set<int>> pt2ins;
    for (int i = 0; i < tetras_num; ++i) {
        // if (_showlog) ofs1 << i << ", tt: ";
        Fixed_alpha_shape_3::Cell_handle tt = interior_tetras[i];
        for (int j = 0; j < 4; ++j) {
            size_t vid = tt->vertex(j)->info();
            v2ins[vid].insert(i);
            //      Point p = tt->vertex(j)->point();
            //      pt2ins[p].insert(i);
            // if (_showlog) ofs1 << p << ", ";
        }
        // if (_showlog) ofs1 << endl;
    } // end tetras
    for (int i = 0; i < tris_num; ++i) {
        // if (_showlog) ofs1 << i + tetras_num << ", tri: ";
        Fixed_alpha_shape_3::Facet tri = singular_tris[i];
        auto tmp_tetra = tri.first;
        for (int j = 0; j < 4; j++) {
            if (j != tri.second) {
                size_t vid = tmp_tetra->vertex(j)->info();
                v2ins[vid].insert(i + tetras_num);
                //        Point p = tmp_tetra->vertex(j)->point();
                //        pt2ins[p].insert(i + tetras_num);
                //        if (_showlog) ofs1 << p << ", ";
            }
        }
        // if (_showlog) ofs1 << endl;
    } // end tris
    for (int i = 0; i < edges_num; ++i) {
        // if (_showlog) ofs1 << i + tetras_num + tris_num << ", edge: ";
        Fixed_alpha_shape_3::Edge e = singular_edges[i];
        auto tmp_tetra = e.get<0>();
        int pid1, pid2;
        pid1 = e.get<1>();
        pid2 = e.get<2>();
        size_t vid1, vid2;
        vid1 = tmp_tetra->vertex(pid1)->info();
        vid2 = tmp_tetra->vertex(pid2)->info();
        v2ins[vid1].insert(i + tetras_num + tris_num);
        v2ins[vid2].insert(i + tetras_num + tris_num);
        //    Point p1, p2;
        //    p1 = tmp_tetra->vertex(pid1)->point();
        //    p2 = tmp_tetra->vertex(pid2)->point();
        //    pt2ins[p1].insert(i + tetras_num + tris_num);
        //    pt2ins[p2].insert(i + tetras_num + tris_num);
        //    if (_showlog) ofs1 << p1 << ", " << p2 << endl;
    } // end edges
    // if (_showlog) ofs1.close();

    cout << "Summary of entire alpha shape\n";
    cout << "# of tetra: " << tetras_num << endl;
    cout << "# of tris: " << tris_num << endl;
    cout << "# of edges: " << edges_num << endl;
    cout << "# of pts: " << v2ins.size() << endl;
    //  cout << "# of pts: " << pt2ins.size() << endl;
    //    for (const auto& p2i : pt2ins) {
    //      cout << p2i.first << " ";
    //      for (const auto& i : p2i.second) {
    //        cout << i << " ";
    //      }
    //      cout << endl;
    //    }

    end_memory = xx_getMemoryInstantUsage();
    cout << "finish vertex 2 instances\n";
    cout << "used memory:" << end_memory << " MB\n";
    cout << "used memory of vertex2instances: " << end_memory - st_memory << " MB" << endl;

    cout << "start to find connected objects\n";
    // std::map<Point, int> visited_pts;  // default is not visited
    std::map<int, int> visited_vs; // default is not visited
    int sgs_num = 1; // valid label starts from 1
    for (const auto &p2i: v2ins) {
        int vid = p2i.first;
        if (visited_vs.find(vid) != visited_vs.end() && visited_vs[vid] > 0) {
            // visited
            continue;
        }
        // a new connected object starts
        // collect tetras, triangles, edges
        // cout << "\na new object starts from point: " << p << endl;
        // std::vector<Point> to_process_pts = {p};
        std::vector<int> to_process_vs = {vid};
        std::set<int> con_ins;
        while (!to_process_vs.empty()) {
            // Point cur_p = to_process_pts.back();
            int cur_vid = to_process_vs.back();
            //      to_process_pts.pop_back();
            to_process_vs.pop_back();
            if (visited_vs.find(cur_vid) != visited_vs.end() && visited_vs[cur_vid] > 0) {
                continue;
            }
            visited_vs[cur_vid] = 1;
            if (_savelbl) {
                // todo: implement it. use info() not point. change the ptlbls type
                // label point: do a little test, whether point is labeled before
                //                if (ptlbls.find(cur_p) != ptlbls.end()) {
                //                  int old_lbl = ptlbls[cur_p];
                //                  if (old_lbl != sgs_num) {
                //                    cout << "ERROR: Point is labeled before\n";
                //                    cout << "Point: " << cur_p << ", old lbl: " << old_lbl << ", current lbl: " <<
                //                    sgs_num << endl; exit(1);
                //                  }
                //                }
                //                ptlbls[cur_p] = sgs_num;
            }

            // cout << "process point: " << cur_p << endl;
            // cout << "# of its instances: " << pt2ins[cur_p].size() << endl;
            for (const auto &ins: v2ins[cur_vid]) {
                // cout << "instance: " << ins << ": ";
                con_ins.insert(ins);
                // find connected points
                if (ins < tetras_num) {
                    // tetra
                    // cout << "tt: ";
                    Fixed_alpha_shape_3::Cell_handle tt = interior_tetras[ins];
                    for (int j = 0; j < 4; ++j) {
                        // Point np = tt->vertex(j)->point();
                        int np = tt->vertex(j)->info();
                        // cout << np << ", ";
                        if (np != cur_vid) {
                            to_process_vs.push_back(np);
                        }
                    }
                    // cout << endl;
                } else if (ins < (tetras_num + tris_num)) {
                    // triangle
                    // cout << "tri: ";
                    Fixed_alpha_shape_3::Facet tri = singular_tris[ins - tetras_num];
                    auto tmp_tetra = tri.first;
                    for (int j = 0; j < 4; j++) {
                        if (j != tri.second) {
                            //              Point np = tmp_tetra->vertex(j)->point();
                            int np = tmp_tetra->vertex(j)->info();
                            // cout << np << ", ";
                            if (np != cur_vid) {
                                to_process_vs.push_back(np);
                            }
                        }
                    }
                    // cout << endl;
                } else {
                    // edge
                    // cout << "edge: ";
                    Fixed_alpha_shape_3::Edge e = singular_edges[ins - tetras_num - tris_num];
                    auto tmp_tetra = e.get<0>();
                    int pid1, pid2;
                    pid1 = e.get<1>();
                    pid2 = e.get<2>();
                    int p1, p2;
                    p1 = tmp_tetra->vertex(pid1)->info();
                    p2 = tmp_tetra->vertex(pid2)->info();
                    //          Point p1, p2;
                    //          p1 = tmp_tetra->vertex(pid1)->point();
                    //          p2 = tmp_tetra->vertex(pid2)->point();
                    if (p1 != cur_vid) {
                        to_process_vs.push_back(p1);
                    }
                    if (p2 != cur_vid) {
                        to_process_vs.push_back(p2);
                    }
                    // cout << p1 << " " << p2;
                    // cout << endl;
                }
            }
        } // finish DFS

        // get the con_ins vector: instance id
        // get real objects from con_ins
        std::vector<Fixed_alpha_shape_3::Cell_handle> cells; // tetra
        std::vector<Fixed_alpha_shape_3::Facet> facets; // triangles
        std::vector<Fixed_alpha_shape_3::Edge> edges; // edges
        for (const auto &ins: con_ins) {
            if (ins < tetras_num) {
                // tetra
                Fixed_alpha_shape_3::Cell_handle tt = interior_tetras[ins];
                cells.push_back(tt);
            } else if (ins < tetras_num + tris_num) {
                // triangle
                Fixed_alpha_shape_3::Facet tri = singular_tris[ins - tetras_num];
                facets.push_back(tri);
            } else {
                // edge
                Fixed_alpha_shape_3::Edge e = singular_edges[ins - tris_num - tetras_num];
                edges.push_back(e);
            }
        } // fill up cells, facets, edges
        // output
        // output all connected objects, filtering work is done separately
        // but let have removed obviously small objects
        int tmp_total_num = 4 * cells.size() + 3 * facets.size() + 2 * edges.size();
        if (tmp_total_num < 20) {
            continue;
        }
        // info
        std::cout << "\nFind a component\nThe simple connected object has: " << std::endl;
        std::cout << cells.size() << " interior tetrahedrons" << std::endl;
        std::cout << facets.size() << " singular triangles" << std::endl;
        std::cout << edges.size() << " singular edges" << std::endl;
        std::stringstream ss;
        ss << outdir << "/s" << sgs_num << "_c" << cells.size() << "f" << facets.size() << "e" << edges.size()
                << ".off";
        string obj_file = ss.str();
        write_object_off(cells, facets, edges, obj_file);
        sgs_num = sgs_num + 1;

        std::cout << "try to clear vectors to save memory\n";
        st_memory = xx_getMemoryInstantUsage();
        // clear vectors to save memory?
        // https://www.techiedelight.com/delete-vector-free-memory-cpp/
        std::vector<int>().swap(to_process_vs);
        std::vector<Fixed_alpha_shape_3::Cell_handle>().swap(cells);
        std::vector<Fixed_alpha_shape_3::Facet>().swap(facets);
        std::vector<Fixed_alpha_shape_3::Edge>().swap(edges);
        std::set<int>().swap(con_ins);
        std::cout << "clean temperature contents: vector, set, etc\n";
        end_memory = xx_getMemoryInstantUsage();
        cout << "memory released: " << end_memory - st_memory << " MB\n";

        // update sgs file
        string sgs_file = outdir + "/sgs_list.txt";
        cout << "update sgs file: " << sgs_file << endl;
        // update sgs list with the absolute path
        ofstream ofs2(sgs_file, std::fstream::app);
        string ab_filepath;
        char *real_path = realpath(obj_file.c_str(), NULL);
        ab_filepath = real_path;
        free(real_path);
        cout << obj_file << " -> " << ab_filepath << endl;
        ofs2 << ab_filepath << endl;
        ofs2.close();
    } // finish visit all points
    cout << "# of components: " << sgs_num - 1 << endl;
} // end of get_connected_objects_3


// start
void generate_fixed_alpha_shape_only(const std::string &ptsfile, const double &alpha, const std::string &offfile) {
    Delaunay_Fixed dt;

    if (ptsfile.find(".pts") != std::string::npos) {
        std::cout << ".pts format\n";
        std::ifstream ifs(ptsfile);
        std::string line;
        while (std::getline(ifs, line)) {
            std::vector<double> coord;
            str2vector<double>(line, coord);
            if (coord.size() >= 3) {
                Point p(coord[0], coord[1], coord[2]);
                dt.insert(p);
            }
        }
        ifs.close();
    } else if (ptsfile.find(".ply") != std::string::npos) {
        std::cout << ".ply format\n";
        happly::PLYData plyIn(ptsfile);

        // Get vertex positions
        std::vector<std::array<double, 3> > vertices = plyIn.getVertexPositions();

        for (const auto &v: vertices) {
            Point p(v[0], v[1], v[2]);
            dt.insert(p);
        }
    } else if (ptsfile.find(".off") != std::string::npos) {
        std::cout << ".off format not implemented yet\n";
        return;
    } else {
        std::cout << "Unsupported format\n";
        return;
    }

    std::cout << "Delaunay is built\n";

    // Build alpha shape
    Fixed_alpha_shape_3 as(dt, alpha);
    std::cout << "Alpha shape is built\n";

    std::cout << "Writing AS3D: " << offfile << std::endl;
    if (ptsfile.find(".pts") != std::string::npos) {
        std::map<Point, int> ptlbls;
        xx_write_fixed_as_off(as, ptlbls, offfile);
    }
    if (ptsfile.find(".ply") != std::string::npos) {
        xx_write_fixed_as_ply_binary(as, offfile);
    }
}

void generate_fixed_alpha_shape(const std::string &ptsfile, const double &alpha, const std::string &offfile,
                                const int &version = -1, string &output_asfile = (string &) "") {
    // version: version of find components in the alpha shape
    // version: 3: use more memory, faster?. default.
    // version: 4: use fewer memory.

    // todo: clean code
    //  change: output file names: comp_list.
    //  add: control output: entire alpha shape, alpha shape components
    //  study/learn how to design a function with reference string and with a default value.
    //      current solution seems the good.

    float st_memory, end_memory;
    st_memory = xx_getMemoryInstantUsage();
    cout << "\nstart generating alpha shape...\n";
    // cout << "used memory: " << st_memory << " MB" << endl;

    xx_Timer mytimer1;
    mytimer1.start();
    // read points to Delaunay
    Delaunay_Fixed dt;
    ifstream ifs(ptsfile);
    string line;
    while (getline(ifs, line)) {
        std::vector<double> coord;
        str2vector<double>(line, coord);
        Point p(coord[0], coord[1], coord[2]);
        dt.insert(p);
    }
    ifs.close();
    cout << "\nDelaunay is built\n";
    end_memory = xx_getMemoryInstantUsage();
    // cout << "used memory: " << end_memory << " MB" << endl;
    cout << "used memory of Delaunay Triangulation: " << (end_memory - st_memory) / 1024.0 << " GB\n";

    st_memory = xx_getMemoryInstantUsage();
    // cout << "used memory: " << st_memory << " MB" << endl;
    //   generate alpha shape
    Fixed_alpha_shape_3 as(dt, alpha);
    //  https://doc.cgal.org/latest/Alpha_shapes_3/classCGAL_1_1Fixed__alpha__shape__3.html#affb181a32b5dfff420a8c33882e276a0
    //  This operation swaps *this and dt, that is dt is an empty triangulation once the fixed alpha shape is built.
    mytimer1.stop();
    std::cout << "\nAlpha shape computed\n";
    cout << "dimension: " << as.dimension() << ", a FIXED alpha value: " << alpha << "\n";
    mytimer1.stop();
    std::cout << "used time: " << mytimer1.getElapsedTimeInSec() << " s\n";
    end_memory = xx_getMemoryInstantUsage();
    // cout << "used memory: " << end_memory << " MB" << endl;
    cout << "used memory of Alpha shape " << (end_memory - st_memory) / 1024.0 << " GB\n\n";

    // output
    std::map<Point, int> ptlbls;
    if (!offfile.empty()) {
        std::cout << "write: " << offfile << std::endl;
        xx_write_fixed_as_off(as, ptlbls, offfile);
        if (output_asfile.empty()) {
            return; // no need to find components
        }
    }

    // find components in the alpha shape
    std::cout << "\nFinding components in the alpha shape\n";
    std::size_t found = ptsfile.find_last_of("/\\");
    std::size_t found2 = ptsfile.find_last_of('.');
    // use the absolute path
    string outdir = ptsfile.substr(0, found2); //
    // string outdir = ptsfile.substr(found + 1, found2 - found - 1);  //
    stringstream ss;
    ss << fixed << setprecision(3) << outdir << "_a" << alpha << "_sgs";
    recreateFolder(ss.str());
    //mkdir(ss.str().c_str(), 0777);
    cout << "create outdir: " << ss.str() << endl;
    // init list of connected object files
    ofstream ofs(ss.str() + "/sgs_list.txt"); //, std::fstream::trunc
    ofs << "";
    ofs.close();

    output_asfile = ss.str() + "/sgs_list.txt";

    // find and generate connected objects
    bool savelbl = false;
    bool showlog = false;
    if (!offfile.empty()) {
        savelbl = true;
    }

    xx_Timer mytimer;

    mytimer.start();
    st_memory = xx_getMemoryInstantUsage();
    // cout << "used memory: " << st_memory << " MB" << endl;
    //  in version 3 & 4, not able to update ptlbls
    get_connected_objects_3(as, ptlbls, ss.str(), savelbl, showlog);
    // cout << "used memory: " << end_memory << " MB" << endl;
    end_memory = xx_getMemoryInstantUsage();
    cout << "memory of getting connected objects " << (end_memory - st_memory) / 1024.0 << " GB\n";
    mytimer.stop();
    cout << "time usage: " << mytimer.getElapsedTimeInSec() << " s.\n";
    cout << "find simple-connected components \n";


    // may not need to do so...
    cout << "\nclear alpha shape\n";
    st_memory = xx_getMemoryInstantUsage();
    cout << "used memory: " << st_memory << " MB" << endl;
    as.clear();
    end_memory = xx_getMemoryInstantUsage();
    cout << "used memory: " << end_memory << " MB" << endl;
    cout << "released memory: " << st_memory - end_memory << " MB\n";
}


#endif  // CGAL_AS_CGAL_FIXED_ALPHASHAPE_H
