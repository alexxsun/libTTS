#pragma once

// --- Standard Library Includes ---
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <queue>

// --- CGAL Includes ---
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>

// --- Local Dependency for PLY ---
#include "../io_support/happly.h" // Assumes happly is available

// An anonymous namespace is used to keep helper functions "private" to this file,
// preventing linker errors if this header is included in multiple places.
namespace {
    // --- CGAL Type Definitions ---
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point = K::Point_3;
    using Vb_Tri = CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K>;
    using Vb_Fixed = CGAL::Fixed_alpha_shape_vertex_base_3<K, Vb_Tri>;
    using Cb_Fixed = CGAL::Fixed_alpha_shape_cell_base_3<K>;
    using Tds_Fixed = CGAL::Triangulation_data_structure_3<Vb_Fixed, Cb_Fixed>;
    using Delaunay_Fixed = CGAL::Delaunay_triangulation_3<K, Tds_Fixed, CGAL::Fast_location>;
    using Fixed_alpha_shape = CGAL::Fixed_alpha_shape_3<Delaunay_Fixed>;

    // --- Internal Helper Functions ---

    /**
     * @brief Loads points from a file (.pts, .txt, .ply) into a vector.
     * @param ptsfile Path to the input point cloud file.
     * @param points Vector to be filled with loaded points.
     * @return True on success, false on failure.
     */
    bool _load_points(const std::filesystem::path &ptsfile, std::vector<Point> &points) {
        if (!std::filesystem::exists(ptsfile)) {
            std::cerr << "Error: Input file not found: " << ptsfile << std::endl;
            return false;
        }

        std::string ext = ptsfile.extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        if (ext == ".pts" || ext == ".txt") {
            std::ifstream ifs(ptsfile);
            if (!ifs.is_open()) return false;
            double x, y, z;
            while (ifs >> x >> y >> z) {
                points.emplace_back(x, y, z);
                ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        } else if (ext == ".ply") {
            try {
                happly::PLYData plyIn(ptsfile.string());
                auto vertices = plyIn.getVertexPositions();
                for (const auto &v: vertices) {
                    points.emplace_back(v[0], v[1], v[2]);
                }
            } catch (const std::exception &e) {
                std::cerr << "Error reading PLY file: " << e.what() << std::endl;
                return false;
            }
        } else {
            std::cerr << "Error: Unsupported input file format: " << ext << std::endl;
            return false;
        }
        return true;
    }

    /**
     * @brief Writes a mesh component to a .off file.
     * @param vertices A vector of the points in the mesh component.
     * @param cells A vector of the tetrahedra in the mesh component.
     * @param facets A vector of the triangles in the mesh component.
     * @param edges A vector of the edges in the mesh component.
     * @param outfile The path for the output .off file.
     * @return True on success, false on failure.
     */
    bool _write_component_off(
        const std::vector<Point> &vertices,
        const std::vector<std::array<size_t, 4> > &cells,
        const std::vector<std::array<size_t, 3> > &facets,
        const std::vector<std::array<size_t, 2> > &edges,
        const std::filesystem::path &outfile) {
        std::ofstream ofs(outfile);
        if (!ofs.is_open()) {
            std::cerr << "Error: Could not open file for writing: " << outfile << std::endl;
            return false;
        }
        // print summary
        std::cout << "vertices #:" << vertices.size() << std::endl;
        std::cout << "tetra #: " << cells.size() << std::endl;
        std::cout << "triangles #: " << facets.size() << std::endl;
        std::cout << "edges #: " << edges.size() << std::endl;

        ofs << "OFF\n";
        ofs << vertices.size() << " " << (cells.size() + facets.size() + edges.size()) << " 0\n";

        ofs << std::fixed << std::setprecision(6);
        for (const auto &p: vertices) {
            ofs << p.x() << " " << p.y() << " " << p.z() << "\n";
        }

        for (const auto &cell: cells) {
            ofs << "4 " << cell[0] << " " << cell[1] << " " << cell[2] << " " << cell[3] << "\n";
        }
        for (const auto &facet: facets) {
            ofs << "3 " << facet[0] << " " << facet[1] << " " << facet[2] << "\n";
        }
        for (const auto &edge: edges) {
            ofs << "2 " << edge[0] << " " << edge[1] << "\n";
        }

        std::cout << "Successfully wrote mesh to:\n" << outfile.string() << std::endl;
        return true;
    }

    /**
     * @brief Writes a mesh component to a binary .ply file.
     * @param vertices A vector of the points in the mesh component.
     * @param cells A vector of the tetrahedra in the mesh component.
     * @param facets A vector of the triangles in the mesh component.
     * @param edges A vector of the edges in the mesh component.
     * @param outfile The path for the output .ply file.
     * @return True on success, false on failure.
     */
    bool _write_component_ply(
        const std::vector<Point> &vertices,
        const std::vector<std::array<size_t, 4> > &cells,
        const std::vector<std::array<size_t, 3> > &facets,
        const std::vector<std::array<size_t, 2> > &edges,
        const std::filesystem::path &outfile) {
        try {
            std::vector<std::array<double, 3> > vertex_positions;
            vertex_positions.reserve(vertices.size());
            for (const auto &p: vertices) {
                vertex_positions.push_back({p.x(), p.y(), p.z()});
            }

            // Combine all simplex types into a single list of faces for PLY
            std::vector<std::vector<size_t> > ply_total_faces;
            ply_total_faces.reserve(cells.size() + facets.size() + edges.size());

            for (const auto &cell: cells) {
                ply_total_faces.emplace_back(cell.begin(), cell.end());
            }
            for (const auto &facet: facets) {
                ply_total_faces.emplace_back(facet.begin(), facet.end());
            }
            for (const auto &edge: edges) {
                ply_total_faces.emplace_back(edge.begin(), edge.end());
            }

            happly::PLYData plyOut;
            plyOut.addVertexPositions(vertex_positions);
            plyOut.addFaceIndices(ply_total_faces);
            plyOut.write(outfile.string(), happly::DataFormat::Binary);

            std::cout << "Successfully wrote " << ply_total_faces.size() << " faces to\n" << outfile.string() <<
                    std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error writing PLY file: " << e.what() << std::endl;
            return false;
        }
        return true;
    }

    /**
     * @brief Finds and saves the connected components of the alpha shape.
     * @param as The computed alpha shape.
     * @param out_dir The directory to save component files.
     * @param format The desired output format for each component.
     */
    // todo: check the function
    void _find_and_save_components(const Fixed_alpha_shape &as, const std::filesystem::path &out_dir,
                                   const std::string &format) {
        std::vector<Fixed_alpha_shape::Cell_handle> interior_cells;
        as.get_alpha_shape_cells(std::back_inserter(interior_cells), Fixed_alpha_shape::INTERIOR);

        std::vector<Fixed_alpha_shape::Facet> singular_facets;
        as.get_alpha_shape_facets(std::back_inserter(singular_facets), Fixed_alpha_shape::SINGULAR);

        std::vector<Fixed_alpha_shape::Edge> singular_edges;
        as.get_alpha_shape_edges(std::back_inserter(singular_edges), Fixed_alpha_shape::SINGULAR);

        // Assign a unique, contiguous ID to each vertex using its built-in info() field.
        size_t vertex_count = 0;
        for (auto vit = as.finite_vertices_begin(); vit != as.finite_vertices_end(); ++vit) {
            vit->info() = vertex_count++;
        }

        // Build adjacency list: vertex ID -> list of adjacent vertex IDs
        std::vector<std::set<size_t> > adj(vertex_count);
        auto add_edge = [&](Fixed_alpha_shape::Vertex_handle v1, Fixed_alpha_shape::Vertex_handle v2) {
            adj[v1->info()].insert(v2->info());
            adj[v2->info()].insert(v1->info());
        };

        for (const auto &cell: interior_cells) {
            for (int i = 0; i < 4; ++i)
                for (int j = i + 1; j < 4; ++j) {
                    add_edge(cell->vertex(i), cell->vertex(j));
                }
        }
        for (const auto &facet: singular_facets) {
            auto v1 = facet.first->vertex((facet.second + 1) % 4);
            auto v2 = facet.first->vertex((facet.second + 2) % 4);
            auto v3 = facet.first->vertex((facet.second + 3) % 4);
            add_edge(v1, v2);
            add_edge(v2, v3);
            add_edge(v3, v1);
        }
        for (const auto &edge: singular_edges) {
            add_edge(edge.first->vertex(edge.second), edge.first->vertex(edge.third));
        }

        // Find connected components using BFS
        std::vector<bool> visited(vertex_count, false);
        int component_count = 0;
        for (size_t i = 0; i < vertex_count; ++i) {
            if (!visited[i]) {
                component_count++;
                std::vector<size_t> current_component_v_ids;
                std::queue<size_t> q;

                q.push(i);
                visited[i] = true;
                while (!q.empty()) {
                    size_t u = q.front();
                    q.pop();
                    current_component_v_ids.push_back(u);
                    for (size_t v: adj[u]) {
                        if (!visited[v]) {
                            visited[v] = true;
                            q.push(v);
                        }
                    }
                }

                // --- Reconstruct the component mesh for saving ---
                std::set<size_t> component_v_set(current_component_v_ids.begin(), current_component_v_ids.end());
                std::map<size_t, size_t> old_id_to_new_id;
                std::vector<Point> component_vertices;
                component_vertices.reserve(current_component_v_ids.size());

                // Create a map from the old global ID to the new local ID for this component
                for (size_t old_id: current_component_v_ids) {
                    old_id_to_new_id[old_id] = component_vertices.size();
                }

                // Build the vertex list for this component
                std::vector<Fixed_alpha_shape::Vertex_handle> id_to_vh(vertex_count);
                for (auto vit = as.finite_vertices_begin(); vit != as.finite_vertices_end(); ++vit) {
                    id_to_vh[vit->info()] = vit;
                }
                for (size_t old_id: current_component_v_ids) {
                    component_vertices.push_back(id_to_vh[old_id]->point());
                }

                std::vector<std::array<size_t, 4> > component_cells;
                std::vector<std::array<size_t, 3> > component_facets;
                std::vector<std::array<size_t, 2> > component_edges;

                // Filter simplices belonging to the current component and map to new local IDs
                for (const auto &cell: interior_cells) {
                    if (component_v_set.count(cell->vertex(0)->info())) {
                        component_cells.push_back({
                            old_id_to_new_id.at(cell->vertex(0)->info()), old_id_to_new_id.at(cell->vertex(1)->info()),
                            old_id_to_new_id.at(cell->vertex(2)->info()), old_id_to_new_id.at(cell->vertex(3)->info())
                        });
                    }
                }
                for (const auto &facet: singular_facets) {
                    if (component_v_set.count(facet.first->vertex((facet.second + 1) % 4)->info())) {
                        component_facets.push_back({
                            old_id_to_new_id.at(facet.first->vertex((facet.second + 1) % 4)->info()),
                            old_id_to_new_id.at(facet.first->vertex((facet.second + 2) % 4)->info()),
                            old_id_to_new_id.at(facet.first->vertex((facet.second + 3) % 4)->info())
                        });
                    }
                }
                for (const auto &edge: singular_edges) {
                    if (component_v_set.count(edge.first->vertex(edge.second)->info())) {
                        component_edges.push_back({
                            old_id_to_new_id.at(edge.first->vertex(edge.second)->info()),
                            old_id_to_new_id.at(edge.first->vertex(edge.third)->info())
                        });
                    }
                }

                // Save the component
                try {
                    if (!std::filesystem::exists(out_dir)) {
                        std::filesystem::create_directories(out_dir);
                    }
                    std::filesystem::path component_path =
                            out_dir / ("component_" + std::to_string(component_count) + "." + format);

                    if (format == "off") {
                        _write_component_off(component_vertices, component_cells, component_facets, component_edges,
                                             component_path);
                    } else if (format == "ply") {
                        _write_component_ply(component_vertices, component_cells, component_facets, component_edges,
                                             component_path);
                    }
                } catch (const std::filesystem::filesystem_error &e) {
                    std::cerr << "Filesystem error: " << e.what() << std::endl;
                }
            }
        }
        std::cout << "Found and processed " << component_count << " connected components." << std::endl;
    }

    /**
     * @brief Dispatches to the correct mesh writing function based on format.
     * @param as The computed alpha shape.
     * @param outfile The path for the output file.
     * @param format The desired format ("off" or "ply").
     * @return True on success, false on failure.
     */
    bool _write_mesh(const Fixed_alpha_shape &as, const std::filesystem::path &outfile, const std::string &format) {
        // This function now acts as a dispatcher for writing the *entire* alpha shape.
        // It gathers all vertices and faces and calls the appropriate component writer.

        // Assign IDs directly using the info() field.
        std::vector<Point> all_vertices;
        all_vertices.reserve(as.number_of_vertices());
        size_t vertex_count = 0;
        for (auto vit = as.finite_vertices_begin(); vit != as.finite_vertices_end(); ++vit) {
            vit->info() = vertex_count++;
            all_vertices.push_back(vit->point());
        }

        std::vector<std::array<size_t, 4> > all_cells;
        std::vector<Fixed_alpha_shape::Cell_handle> interior_cells;
        as.get_alpha_shape_cells(std::back_inserter(interior_cells), Fixed_alpha_shape::INTERIOR);
        for (const auto &cell: interior_cells) {
            all_cells.push_back({
                cell->vertex(0)->info(), cell->vertex(1)->info(), cell->vertex(2)->info(), cell->vertex(3)->info()
            });
        }

        std::vector<std::array<size_t, 3> > all_facets;
        /*
         * https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#ad6a20b45e66dfb690bfcdb8438e9fcae
         * (c,i) is the facet of c opposite to the vertex of index i.
         * c: cell <- tetrahedron; i: vertex
         */
        std::vector<Fixed_alpha_shape::Facet> singular_facets;
        as.get_alpha_shape_facets(std::back_inserter(singular_facets), Fixed_alpha_shape::SINGULAR);
        for (const auto &facet: singular_facets) {
            all_facets.push_back({
                facet.first->vertex((facet.second + 1) % 4)->info(),
                facet.first->vertex((facet.second + 2) % 4)->info(),
                facet.first->vertex((facet.second + 3) % 4)->info()
            });
        }

        std::vector<std::array<size_t, 2> > all_edges;
        /*
         * https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#af31db7673a6d7d28c0bb90a3115ac695
         * (c,i,j) is the edge of cell c whose vertices indices are i and j.
         */
        std::vector<Fixed_alpha_shape::Edge> singular_edges;
        as.get_alpha_shape_edges(std::back_inserter(singular_edges), Fixed_alpha_shape::SINGULAR);
        for (const auto &edge: singular_edges) {
            all_edges.push_back({edge.first->vertex(edge.second)->info(), edge.first->vertex(edge.third)->info()});
        }

        if (format == "off") {
            return _write_component_off(all_vertices, all_cells, all_facets, all_edges, outfile);
        }
        if (format == "ply") {
            return _write_component_ply(all_vertices, all_cells, all_facets, all_edges, outfile);
        }
        std::cerr << "Error: Unsupported output format '" << format << "'" << std::endl;
        return false;
    }
} // end anonymous namespace

// --- Public API Functions ---

/**
 * @brief Generates a 3D alpha shape and saves it to a single file.
 * @param ptsfile_str Path to the input point cloud file.
 * @param alpha_sq The squared alpha value for the shape computation.
 * @param outfile_str Path for the output mesh file.
 * @param out_format The desired output format ("off" or "ply").
 * @return 1 on success, 0 on failure.
 */
int generate_fixed_alpha_shape_only_v2(
    const std::string &ptsfile_str,
    const double &alpha_sq,
    const std::string &outfile_str,
    const std::string &out_format) {
    std::vector<Point> points;
    if (!_load_points(ptsfile_str, points)) {
        return 0; // Failure
    }

    if (points.empty()) {
        std::cerr << "Error: No points loaded from file." << std::endl;
        return 0;
    }

    Delaunay_Fixed dt(points.begin(), points.end());
    Fixed_alpha_shape as(dt, alpha_sq);

    if (!_write_mesh(as, outfile_str, out_format)) {
        return 0; // Failure
    }

    return 1; // Success
}

/**
 * @brief Generates a 3D alpha shape and extracts its connected components.
 * @param ptsfile_str Path to the input point cloud file.
 * @param alpha_sq The squared alpha value for the shape computation.
 * @param out_asfile Path to save the complete alpha shape mesh. Can be empty.
 * @param out_component_dir Path to the directory where component files will be saved.
 * @param out_format The desired output format for all files ("off" or "ply").
 * @return 1 on success, 0 on failure.
 */
int generate_fixed_alpha_shape_v2(
    const std::string &ptsfile_str,
    const double &alpha_sq,
    std::string &out_asfile,
    std::string &out_component_dir,
    const std::string &out_format) {
    std::vector<Point> points;
    if (!_load_points(ptsfile_str, points)) {
        return 0; // Failure
    }

    if (points.empty()) {
        std::cerr << "Error: No points loaded from file." << std::endl;
        return 0;
    }

    Delaunay_Fixed dt(points.begin(), points.end());
    Fixed_alpha_shape as(dt, alpha_sq);

    // Save the main alpha shape if a path is provided
    if (!out_asfile.empty()) {
        if (!_write_mesh(as, out_asfile, out_format)) {
            std::cerr << "Error: Failed to save the main alpha shape." << std::endl;
            // Continue to component extraction even if this fails, or return 0
        }
    }

    // Find and save components
    if (!out_component_dir.empty()) {
        _find_and_save_components(as, out_component_dir, out_format);
    }

    return 1; // Success
}
