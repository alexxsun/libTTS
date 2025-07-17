#ifndef SIMPLICIALCOMPLEX_H
#define SIMPLICIALCOMPLEX_H

// https://stackoverflow.com/questions/53203970/why-boostbind-insists-pulling-boostplaceholders-into-global-namespace
//#include <boost/bind.hpp>
#include <boost/bind/bind.hpp>
#include <boost/function.hpp>
#include <cmath>
#include <forward_list>
#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "topsimplex.h"
#include "vertex.h"

class SimplicialComplex {
protected:
    // simplex dim: dim id set by our code
    map<int, int> realIndex;
    // dim id: top simplexes. Note: the dim id is not dim of simplex. see realIndex[dim] = dim_id
    vector<vector<TopSimplex> > topSimplexes;

    vector<forward_list<int> > adjRelations;

    // Vertex: x y z f1 f2 ...
    vector<Vertex> vertices;

    // empty by default. Call storeFullStart() to encode the VTop relations and boost computation speed.
    vector<vector<explicitS> > topPerVertex;

    // EmbeddingSpace == number of coordinates for each vertex
    // SimplicialComplexSpace == dimension of the biggest top simplex

public:
    SimplicialComplex();

    // points extent.
    double sc_min_x, sc_min_y, sc_min_z, sc_max_z;

    size_t getVertexCoordSize();

    //------------------------------------------------------------------------------------------
    // functions for reading input files. Supported formats [OFF,GMV]
    void readOFF(const char *file); // read .OFF file
    void readPLY(const char *file); // read .PLY file
    void readOFF(const char *file, const int &funID);

    void readTS(const char *file); // read .TS  file
    void readGMV(const char *file); // read .GMV file
    void readSimplexes(char *file); // TODO //set of simplexes. Vertices have no coordinates
    void readPoints(char *name_in_file, float threshold,
                    uint dimension); // set of points on which a Vietoris-Rips
    // complex is computed
    void readGraph(char *name_in_file); // set of arcs of a graph based on which clicks are computed
    void readIA(char *file); // read IA file format
    void readSquareGrid(const char *file, int xRes, int yRes);

    void readCubicalGrid(const char *file, int xRes, int yRes, int zRes);

    // writing files
    void saveIA(char *file); // write IA file
    void saveVTK(char *file); // write VTK file


    // core functions

    // explicitly store/remove the complete set of top-simplexes for each vertex
    // (speedup computation bur requires more memory)
    void storeFullStar();

    void emptyFullStar();

    // return: the number of vertices
    int getVerticesNum();

    // return: the number of top-simplexes of dimension d
    int getTopSimplexesNum(int d);

    // return: the simplicial complex dimension (dimension of its maximal simplex)
    int getComplexDim();

    // return: return the number of top-simplexes of each dimension (dimensions having 0 top-simplexes are ignored)
    vector<int> getTopSimplexesSet();

    // return: vertex of index i
    Vertex &getVertex(int i);

    // return: set of top-simplexes of dimendion d
    vector<TopSimplex> &getTopSimplices(int d);

    // return: top-simplex corresponding to the input simplex encoded in the explicit way
    TopSimplex &getTopSimplex(const explicitS);


    //------------------------------------------------------------------------------------------
    // functions for simplexes explicitly encoded in the IA* data structure
    // (explicitS)

    // return: implicit representation of a simplex computed from its explicit representation
    implicitS toImplicit(const explicitS &);

    // return: the set of simplexes in the link of s.
    // NOTE: only the simplexes of highets dimension are returned, i.e. if edge e
    // is in the link of s, e is returned by the function but its vertices are not.
    vector<implicitS> *link(const explicitS &);

    // return: simplexes of dimension k on the boundary of s
    vector<implicitS> *boundaryk(const explicitS &s, uint k);

    // return: simplexes of dimension k on the coboundary of s
    vector<implicitS> *coboundaryk(const explicitS &s, uint k);

    // return: simplexes adjacent to s and sharing the i-th face of s
    vector<explicitS> *topAdjacent(const explicitS &s, uint i);

    // return: top-simplexes on the coboundary of s (or, in other words, incident to s)
    // NOTE: valid only for vertices (a top-simplex cannot have simplexes on their coboundary)
    vector<explicitS> *topStar(const explicitS &s);

    // return: top-simplexes of dimension d on the coboundary of s (or, in other
    // words, incident to s) NOTE: valid only for vertices (top-simplex cannot have
    // simplexes on their coboundary)
    vector<explicitS> *topStar(const explicitS &s, int d);

    //------------------------------------------------------------------------------------------
    // functions for simplexes not encoded in the IA* and implicitly represented
    // (implicitS) [see implicitS.h for details]


    // return: explicit representation of a simplex computed from its implicit
    // representation NOTE: only a vertex or a top simplex have a valid explicit
    // representation;
    explicitS toExplicit(const implicitS &);

    // return: true if s1 and s2 represent the same simplex
    // NOTE: use theSame() when you do not know if s1 and s2 are represented based
    // on the same top-simplex. Otherwise operator== is faster. [see implicitS.h
    // for details]
    bool theSame(const implicitS &s1, const implicitS &s2);

    // return: the set of simplexes in the link of s.
    // NOTE: only the simplexes of highest dimension are returned, i.e. if edge e
    // is in the link of s, e is returned by the function but its vertices are not.
    vector<implicitS> *link(const implicitS &s);

    // return: simplexes of dimension k on the boundary of s
    vector<implicitS> *boundaryk(const implicitS &s, uint k);

    // return: simplexes of dimension k on the coboundary of s
    vector<implicitS> *coboundaryk(const implicitS &, uint);

    // return: simplexes adjacent to s
    vector<implicitS> *adjacents(const implicitS &s);

    // return: top-simplexes on the coboundary of s (or, in other words, incident to s)
    vector<explicitS> *topStar(const implicitS &s);

    // return: true if "small" is on the boundary of "big"
    bool areIncident(const implicitS &small, const implicitS &big);

    //------------------------------------------------------------------------------------------
    // functions for geometry, vectors, etc
    Vertex barycenter(const implicitS &s);

    Vertex barycenter(const explicitS &s);

    vector<float> vectorSum(const vector<float> &v1, const vector<float> &v2);

    vector<float> vectorSubtr(const vector<float> &v1, const vector<float> &v2);

    vector<float> scalarProd(const vector<float> &v1, const vector<float> &v2);

    vector<float> prod(const vector<float> &v1, float prod);

    float dotProd(const vector<float> &v1, const vector<float> &v2);

    float dist(const vector<float> &v1, const vector<float> &v2);

    float norm(const vector<float> &v1);

protected:
    void buildDataStructure();

    // function for initializing the IA*
    void buildDataStructure_parallel();

    // return: implicit representation of a vertex (inside top simplex simpl)
    // computed from its explicit representation
    implicitS toImplicitInsideTop(const explicitS vertex, const explicitS simpl);

    vector<implicitS> *inTop(const explicitS &topSimpl, uint dim);

    uint32_t nextSubset(uint32_t v);

    vector<int> filter(const vector<int> &v, uint32_t mask);

    vector<implicitS> *getSubsets(const vector<int> &arr, uint32_t k);

    //(given a vertex v, and a top simplex of dimension k incident in v, retrieve
    // all the k-1 connected simplexes still incident in v)
    forward_list<explicitS> *incidentCluster(const explicitS v, const explicitS s);

    // recursive function used for computing boundary and coboundary simplexes.
    void recursive_insert(implicitS simplex, const implicitS &original, uint pos, uint dim,
                          forward_list<implicitS> *ret);


    // Function for computing the top simplexes given a graph
    void build_top_simplex(set<uint> setR, set<uint> setP, set<uint> setX, const vector<set<uint> *> &arcs,
                           map<int, list<TopSimplex> *> *top_simplexes_local);
};

#endif  // SIMPLICIALCOMPLEX_H
