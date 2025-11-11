#include "simplicialcomplex.h"
#include <omp.h>
#include <unordered_set>
#include <chrono>


SimplicialComplex::SimplicialComplex() {
    completeCoboundaryTop.clear();
}

int SimplicialComplex::getComplexDim() { return getTopSimplexesSet().back(); }

int SimplicialComplex::getVerticesNum() { return vertices.size(); }

size_t SimplicialComplex::getVertexCoordSize() {
    // fixme: do all vertices have the same diemsions?
    return vertices[0].getCoordinates().size();
}

int SimplicialComplex::getTopSimplexesNum(int dim) {
    if (realIndex.find(dim) == realIndex.end())
        return 0;
    else
        return topSimplexes[realIndex[dim]].size();
}

vector<int> SimplicialComplex::getTopSimplexesSet() {
    vector<int> result;
    for (auto r: realIndex) {
        result.push_back(r.first);
    }
    return result;
}

void SimplicialComplex::buildDataStructure() {
    topPerVertex = vector<vector<explicitS> >();

    //  topSimplexes is vector<vector<TopSimplex> >
    for (uint i = 0; i < topSimplexes.size(); i++) {
        int dim = topSimplexes[i][0].getDimension(); // dim of top simplex

        // Note: Adjacency relations building removed - no longer needed
        // We use completeCoboundaryTop instead for incidentCluster()

        // xx: here, we build partial_coboundaryTop for each vertex
        // possible to do the parallelizing?
        vector<set<int> > incidentTop(vertices.size(), set<int>());
        // vertex_index: set of top simplexes including the vertex
        // todo: use entire vertices.size()? seems costly?
        // build partial incidence relations for vertices
        for (uint j = 0; j < topSimplexes[i].size(); j++) {
            TopSimplex tS = topSimplexes[i][j];
            for (int v = 0; v < tS.getDimension() + 1; v++) {
                incidentTop[tS.getVertexIndex(v)].insert(j);
            }
        }

        // Store complete vertex-to-top-simplex mapping
        // OPTIMIZED: We no longer need to build partialCoboundaryTop since topStar() uses completeCoboundaryTop directly
        auto start_storage = std::chrono::high_resolution_clock::now();
        if (completeCoboundaryTop.size() <= vertices.size()) {
            completeCoboundaryTop.resize(vertices.size());
        }
        for (uint j = 0; j < vertices.size(); j++) {
            if (!incidentTop[j].empty()) {
                completeCoboundaryTop[j][dim].insert(incidentTop[j].begin(), incidentTop[j].end());
            }
        }
        auto end_storage = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> storage_elapsed = end_storage - start_storage;
        std::cout << "      Complete mapping storage took: " << storage_elapsed.count() << " seconds\n";

        // NOTE: Removed partialCoboundaryTop building phase - no longer needed!
        // topStar() now uses completeCoboundaryTop directly, so we don't need to:
        // - Call incidentCluster() to find connected components
        // - Store representatives in partialCoboundaryTop
        // This saves significant time during buildDataStructure()
    }
}

#include <chrono>

void SimplicialComplex::buildDataStructure_parallel() {
    using Clock = std::chrono::high_resolution_clock;
    auto start_total = Clock::now();

    // Clear and pre-allocate data structures
    topPerVertex.clear();

    // // set omp threads num
    // int xx_threads = 4;
    // omp_set_num_threads(xx_threads);

    // Process each dimension of top simplexes (sequential)
    for (uint i = 0; i < topSimplexes.size(); i++) {
        auto start_phase = Clock::now();

        if (topSimplexes[i].empty()) continue;

        const int dim = topSimplexes[i][0].getDimension();
        const size_t simplices_count = topSimplexes[i].size();
        const size_t faces_per_simplex = dim + 1;
        const size_t total_faces = simplices_count * faces_per_simplex;

        cout << "dim: " << dim << ", simplex #: " << simplices_count << ", faces #: " << total_faces << endl;

        // Note: Adjacency relations building removed - no longer needed
        // We use completeCoboundaryTop instead for incidentCluster()

        // Build vertex incidence relations
        auto st_phase = Clock::now();
        vector<unordered_set<int> > incidentTop(vertices.size());

        // OPTIMIZED: Use thread-local accumulation to eliminate critical section bottleneck
        // Each thread processes top simplexes independently and accumulates results locally
        // Then merges once per thread (minimal critical sections)
#pragma omp parallel
        {
            // Each thread has its own local storage - no locking needed during accumulation
            vector<unordered_set<int> > local_incidentTop(vertices.size());
            
            // Process top simplexes in parallel - no critical sections in this loop
#pragma omp for nowait
            for (uint j = 0; j < simplices_count; j++) {
                TopSimplex tS = topSimplexes[i][j];
                for (int v = 0; v < tS.getDimension() + 1; v++) {
                    int vertexIdx = tS.getVertexIndex(v);
                    local_incidentTop[vertexIdx].insert(j);  // No lock needed - thread-local!
                }
            }
            
            // Merge thread-local results into global incidentTop
            // Only one critical section per thread (much less contention)
#pragma omp critical
            {
                for (uint v = 0; v < vertices.size(); v++) {
                    if (!local_incidentTop[v].empty()) {
                        incidentTop[v].insert(local_incidentTop[v].begin(), 
                                             local_incidentTop[v].end());
                    }
                }
            }
        }

        // Store complete vertex-to-top-simplex mapping
        // OPTIMIZED: We no longer need to build partialCoboundaryTop since topStar() uses completeCoboundaryTop directly
        auto start_storage = Clock::now();
        if (completeCoboundaryTop.size() <= vertices.size()) {
            completeCoboundaryTop.resize(vertices.size());
        }
        // Parallelize storage - each vertex is independent
#pragma omp parallel for
        for (uint j = 0; j < vertices.size(); j++) {
            if (!incidentTop[j].empty()) {
                completeCoboundaryTop[j][dim].insert(incidentTop[j].begin(), incidentTop[j].end());
            }
        }
        auto end_storage = Clock::now();
        std::chrono::duration<double> storage_elapsed = end_storage - start_storage;
        std::cout << "      Complete mapping storage took: " << storage_elapsed.count() << " seconds\n";

        // debug: find average number of top simplexes per vertex
        vector<int> top_nums;
        for (uint j = 0; j < vertices.size(); j++) {
            if (incidentTop[j].empty()) continue;
            top_nums.push_back(incidentTop[j].size());
        }
        double num_ratio = 1. * top_nums.size() / vertices.size();
        cout << "      to process vertex#: " << top_nums.size() << "/" << vertices.size() << ": " << num_ratio << endl;

        double avg_top_num = 0;
        for (const auto &it: top_nums) {
            avg_top_num += it;
        }
        avg_top_num /= top_nums.size();
        cout << "      avg top per vertex: " << avg_top_num << endl;

        // NOTE: Removed partialCoboundaryTop building phase - no longer needed!
        // topStar() now uses completeCoboundaryTop directly, so we don't need to:
        // - Call incidentCluster() to find connected components
        // - Store representatives in partialCoboundaryTop
        // This saves significant time during buildDataStructure() - no more expensive cluster expansion!
        
        auto ed_phase = Clock::now();
        std::chrono::duration<double> elapsed1 = ed_phase - st_phase;
        std::cout << "      Phase " << i << " (without partialCoboundaryTop building) took: " << elapsed1.count() <<
                " seconds\n";

        auto end_phase = Clock::now();
        std::chrono::duration<double> elapsed = end_phase - start_phase;
        std::cout << "Phase " << i << " total took: " << elapsed.count() << " seconds\n\n";
    }

    auto end_total = Clock::now();
    std::chrono::duration<double> total = end_total - start_total;
    std::cout << "Total time: " << total.count() << " seconds\n";
}

void SimplicialComplex::build_top_simplex(set<uint> setR, set<uint> setP, set<uint> setX,
                                          const vector<set<uint> *> &arcs,
                                          map<int, list<TopSimplex> *> *top_simplexes_local) {
    if (setP.size() == 0 && setX.size() == 0) {
        // set setR as maximal

        if (setR.size() > 1) {
            vector<int> new_top(setR.begin(), setR.end());
            sort(new_top.begin(), new_top.end());
            TopSimplex top(new_top);

            if (top_simplexes_local->find(top.getDimension()) == top_simplexes_local->end()) {
                (*top_simplexes_local)[top.getDimension()] = new list<TopSimplex>();
            }
            (*top_simplexes_local)[top.getDimension()]->push_back(top);
        }
    } else if (setP.size() > 0) {
        // chose here the pivot vertex for P
        vector<uint> vectorP(setP.begin(), setP.end());

        uint pivot = vectorP[0];
        for (int i = 0; i < vectorP.size(); i++) {
            if (arcs[vectorP[i]]->size() > arcs[pivot]->size()) pivot = vectorP[i];
        }

        for (int i = 0; i < vectorP.size(); i++) {
            uint ui = vectorP[i];
            if (arcs[ui]->find(pivot) == arcs[ui]->end()) {
                setP.erase(ui);
                set<uint> adjui(arcs[ui]->begin(), arcs[ui]->end());

                set<uint> newR = setR;
                newR.insert(ui);

                set<uint> newP;
                set_intersection(setP.begin(), setP.end(), adjui.begin(), adjui.end(),
                                 std::inserter(newP, newP.begin()));

                set<uint> newX;
                set_intersection(setX.begin(), setX.end(), adjui.begin(), adjui.end(),
                                 std::inserter(newX, newX.begin()));

                build_top_simplex(newR, newP, newX, arcs, top_simplexes_local);

                setX.insert(ui);
            }
        }
    }
}

Vertex &SimplicialComplex::getVertex(int vertex) { return vertices[vertex]; }

vector<TopSimplex> &SimplicialComplex::getTopSimplices(int dim) {
    if (realIndex.find(dim) == realIndex.end())
        assert(false);
    else
        return topSimplexes[realIndex[dim]];
}

TopSimplex &SimplicialComplex::getTopSimplex(explicitS simpl) {
    assert(simpl.getDim() > 0);
    if (realIndex.find(simpl.getDim()) == realIndex.end())
        assert(false);
    return topSimplexes[realIndex[simpl.getDim()]][simpl.getIndex()];
}

forward_list<explicitS> *SimplicialComplex::incidentCluster(explicitS vertex, explicitS topS) {
    // Implementation using completeCoboundaryTop (no adjRelations needed)
    // topsimplex must be incident to vertex

    forward_list<explicitS> *ret = new forward_list<explicitS>();
    set<explicitS> visited;
    queue<explicitS> adjacentSimplexes;
    
    adjacentSimplexes.push(topS);
    visited.insert(topS);
    
    while (!adjacentSimplexes.empty()) {
        explicitS current = adjacentSimplexes.front();
        TopSimplex &top = getTopSimplex(current);
        int dim = top.getDimension();
        
        // For each face (excluding faces containing the vertex)
        for (int i = 0; i < top.get_nVertices(); i++) {
            if (top.getVertexIndex(i) == vertex.getIndex()) {
                continue; // Skip faces containing the vertex
            }
            
            // Build face vertices (all vertices except the i-th one)
            set<int> faceVertices;
            for (int j = 0; j < top.get_nVertices(); j++) {
                if (j != i) {
                    faceVertices.insert(top.getVertexIndex(j));
                }
            }
            
            // Find all top simplexes of same dimension incident to vertex that share this face
            // Use completeCoboundaryTop[vertex.getIndex()][dim] to get candidates
            if (vertex.getIndex() >= completeCoboundaryTop.size()) {
                continue; // Safety check
            }
            auto &vertexMap = completeCoboundaryTop[vertex.getIndex()];
            auto dimIt = vertexMap.find(dim);
            if (dimIt == vertexMap.end()) {
                continue; // No top simplexes of this dimension incident to vertex
            }
            
            for (int candidateIdx : dimIt->second) {
                explicitS candidate(dim, candidateIdx);
                
                if (visited.find(candidate) != visited.end()) {
                    continue; // Already processed
                }
                
                TopSimplex &candidateTop = getTopSimplex(candidate);
                vector<int> &candidateVertices = candidateTop.getVertices();
                
                // Check if candidate shares the face (has all face vertices)
                bool sharesFace = true;
                for (int fv : faceVertices) {
                    bool found = false;
                    for (int cv : candidateVertices) {
                        if (cv == fv) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        sharesFace = false;
                        break;
                    }
                }
                
                if (sharesFace) {
                    visited.insert(candidate);
                    adjacentSimplexes.push(candidate);
                }
            }
        }
        
        adjacentSimplexes.pop();
    }
    
    ret->insert_after(ret->before_begin(), visited.begin(), visited.end());
    return ret;
}

vector<explicitS> *SimplicialComplex::topStar(const explicitS &vertex) {
    assert(vertex.getDim() == 0);

    // Check cache first
    if (topPerVertex.size() > vertex.getIndex() && topPerVertex[vertex.getIndex()].size() != 0) {
        return new vector<explicitS>(topPerVertex[vertex.getIndex()]);
    }

    // OPTIMIZED: Direct lookup from completeCoboundaryTop - MUCH FASTER!
    // No need for incidentCluster() calls - we already have all top simplexes stored
    vector<explicitS> ret;
    int vertexIdx = vertex.getIndex();
    
    if (vertexIdx < completeCoboundaryTop.size()) {
        auto &vertexMap = completeCoboundaryTop[vertexIdx];
        for (const auto &dimPair : vertexMap) {
            int dim = dimPair.first;
            for (int topIdx : dimPair.second) {
                ret.push_back(explicitS(dim, topIdx));
            }
        }
    }

    return new vector<explicitS>(ret);
}

vector<explicitS> *SimplicialComplex::topStar(const explicitS &vertex, int dimension) {
    // OPTIMIZED: Direct lookup from completeCoboundaryTop - MUCH FASTER!
    // No need for incidentCluster() calls - we already have all top simplexes stored
    vector<explicitS> ret;
    int vertexIdx = vertex.getIndex();
    
    if (vertexIdx < completeCoboundaryTop.size()) {
        auto &vertexMap = completeCoboundaryTop[vertexIdx];
        auto dimIt = vertexMap.find(dimension);
        if (dimIt != vertexMap.end()) {
            for (int topIdx : dimIt->second) {
                ret.push_back(explicitS(dimension, topIdx));
            }
        }
    }

    return new vector<explicitS>(ret);
}

vector<implicitS> *SimplicialComplex::boundaryk(const explicitS &simplex, uint dim) {
    if (simplex.getDim() < dim) {
        printf("No simplexes of dimensions %d on a %d-simplex", dim, simplex.getDim());
        return new vector<implicitS>();
    }

    vector<implicitS> *sset = getSubsets(getTopSimplex(simplex).getVertices(), dim);

    return sset;
}

vector<implicitS> *SimplicialComplex::coboundaryk(const explicitS &simplex, uint dim) {
    return coboundaryk(toImplicit(simplex), dim);
}

vector<implicitS> *SimplicialComplex::inTop(const explicitS &topSimpl, uint dim) {
    vector<implicitS> *vec = new vector<implicitS>();
    if (dim == topSimpl.getDim()) {
        vec->push_back(implicitS(getTopSimplex(topSimpl).getVertices()));
        return vec;
    } else if (dim > topSimpl.getDim())
        return new vector<implicitS>();
    else {
        cout << "XX: what to do here?\n";
        exit(1);
    }
}

void SimplicialComplex::storeFullStar() {
    topPerVertex = vector<vector<explicitS> >(getVerticesNum());
    
    // OPTIMIZED: Direct copy from completeCoboundaryTop - MUCH FASTER!
    // No need to call topStar() which would call incidentCluster() - direct lookup instead
    // Parallelize vertex processing - each vertex is processed independently
    // Using dynamic scheduling because different vertices may have different
    // numbers of top simplexes (workload imbalance)
#pragma omp parallel for schedule(dynamic) num_threads(6)
    for (int i = 0; i < getVerticesNum(); i++) {
        vector<explicitS> tops;
        
        if (i < completeCoboundaryTop.size()) {
            auto &vertexMap = completeCoboundaryTop[i];
            for (const auto &dimPair : vertexMap) {
                int dim = dimPair.first;
                for (int topIdx : dimPair.second) {
                    tops.push_back(explicitS(dim, topIdx));
                }
            }
        }
        
        sort(tops.begin(), tops.end());
        topPerVertex[i] = tops;
    }
}

void SimplicialComplex::emptyFullStar() { topPerVertex.clear(); }
