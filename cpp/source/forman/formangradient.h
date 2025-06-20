#ifndef FORMANGRADIENT_H
#define FORMANGRADIENT_H

#include <cmath>
#include <iomanip>  // std::setprecision
#include <stack>

#include "gradientencoding.h"
#include "morseincidencegraph.h"

// todo: replace boost::function
typedef set<implicitS, boost::function<bool(const implicitS &, const implicitS &)>> SSet;
typedef map<implicitS, Node *, boost::function<bool(const implicitS &, const implicitS &)>> SMap;
typedef map<implicitS, unsigned, boost::function<bool(const implicitS &, const implicitS &)>> SUMap;
typedef map<implicitS, bool, boost::function<bool(const implicitS &, const implicitS &)>> SBMap;

struct Simpl {
public:
    Arc *arc;
    double val;

    inline Simpl(Arc *a, double v) : arc(a), val(v) {}
};

struct SortSimpl {
    bool operator()(const Simpl &s1, const Simpl &s2) {
        if (s1.val == s2.val) {
            return s1.arc < s2.arc;
        }
        return s1.val > s2.val;
    }
};

typedef priority_queue<Simpl, vector<Simpl>, SortSimpl> SimplQ;

class FormanGradient {
protected:  // private:
    // SimplicialComplex sc need to be accessed in FormanGradient's child class
    //  private ->  protected

    vector<uint> filtration;             // for each vertex its filtration value
    vector<vector<float>> scalarValues;  // for each vertex its field value
    // [Vertices x number of fields] // todo: only one field
    vector<vector<uint>> componentBasedFiltration;  // injective function for each component
    // [number of fields x Vertices ]

    GradientEncoding gradient;
    map<uint, SSet> criticalS;  // dim: simplexes

    SimplicialComplex sc;

public:
    FormanGradient(int, char **);

    FormanGradient(const string &infile);

    FormanGradient(const string &infile, const int &funID);

    ~FormanGradient();

    // compute the Forman gradient on the input dataset, if batchTop=true the
    // relation among vertices and incident top simplices is computed once and
    // stored explicitly (higher memory consumption but lower timings)
    void computeFormanGradient(bool batchTop);


    // Visualization functions
    void visMorse();

    // Xin's visualization
    void xx_vis_CriticalCells_01(const string &vtkfile, const bool &scaled = false);

    void xx_vis_VPath_01(const string &vtkfile, const bool &scaled = false, const bool &outdebugfile = false);
    void xx_vis_VPath_12(const string &vtkfile, const bool &scaled = false, const bool &outdebugfile = false);

    void xx_visMorse(const string &output_prefix);


    void xx_critical_cells_net_2d(const string &vtkfile, const string&txtfile);


    void xx_output_critical_cells(const string &outfile);

    int xx_output_critical_cells_vtk(const string &vtkfile);

    void xx_vis_paired_arrows(const string &outfile);

    // todo: check code. not finished. need debug.
    // https://github.com/IuricichF/FormanGradient2D/blob/master/source/LibForman/io.cpp#L221
    void writeVTK_gradient_2d(const string &vtkfile);

    void write_critical_persistence_2d(const string &outfile);


    // private:
    // some functions need to be accessed in FormanGradient's child class,so change private to protected
protected:
    // xx: get saddle (edge)'s connected minima (vertex)
    void saddle2minima(implicitS const &saddle, vector<implicitS> &minima);
    void saddle2maxima(implicitS const &saddle, vector<implicitS> &maxima);



    // xx: get critical tetra's connected triangles
    void critical_tri2tetras(implicitS const &tri, vector<implicitS> &tetras);

    void critical_tetra2tris(const implicitS &tetra, vector<implicitS> &tris);

    // xx: v-path from saddle to min.
    // two paths. each path last item is the min point id (based on the global point id)
    void saddle2min_vpaths(implicitS const &saddle, std::vector<std::vector<int>> &vpaths);
    void saddle2max_vpaths(implicitS const &saddle, std::vector<std::vector<implicitS>> &vpaths);

    SSet *vertexLowerStar(uint vert, uint d);  // compute the lower star of vert (only simplices of dimension=d)

    void splitVertexLowerStar(int v, vector<SSet> &lwStars);
    // split the lower star of a vertex v according to the sublevelsets of the function

    void homotopy_expansion(SSet &);  // apply homotopy expansion on a set of simplices

    int numPairableLowerStar(const implicitS &next, const SSet &sset, implicitS &pair);
    // return the number of simplices pairable with next in the lower star, pair is one of these

    bool isPaired(const implicitS &simpl);  // true if simpl is paired with another simplex

    void setPair(const implicitS &next, const implicitS &pair);
    // set the new gradient pair between next and pair (NOTE: next has to be bigger than pair)

    bool getPair(const implicitS &simpl, implicitS &next);  // next is the simplex paired with simpl

    void freePair(const implicitS &next, const implicitS &pair);  // remove pair (next,pair) from the gradient
    // (NOTE: next has to be bigger than pair)

    vector<uint> simplexFiltration(const implicitS &simpl);
    // return the vector-valued filtration for a simplex. Each component is obtained as the
    // maximum of the filtrations of its vertices

    vector<float> simplexScalarValue(const implicitS &simpl);  // return the vector-valued function for a
    // simplex. Each component is obtained as the maximum of the function values of its vertices

    // Output functions
    void computeDescendingCell(bool output, implicitS const &cell, SSet &desCells);

    void computeAscendingCell(bool output, implicitS const &cell, SSet &desCells);

    void computeAscendingCell_xx(bool output, implicitS const &cell, SSet &desCells);

    void print_out(const char *fileName, list <SSet> const &cells, int param, int dim);

    void out3cells(const list <SSet> &cells, const string &outfile = "descending3cells.vtk");

    void out2cells(const list <SSet> &cells, bool desc);

    void out1cells(const list <SSet> &cells, bool desc);

    void out0cells(const list <SSet> &cells, const string &filename = "");

    void outCriticalPoints(const SSet &cells);

    void outCriticalPoints_01(const SSet &cells);

    void accurate_asc1cells(const list <SSet> &cells, const string&vtkfile="");



    // Compare the filtration index of the two simplexes (single value filtration)
    /* Each simplex is basically a vector of integers.
     * We want to guarantee that each simplex is stored in the path at most once
     * for maintaining storage cost low. Sorted collection needs O(log n) time
     * complexity to check if a new simplex should be inserted or not. Unsorted
     * collection need O(n) time complexity
     */
    bool cmpSimplexesFiltr(const implicitS &lhs, const implicitS &rhs) {
        if (lhs.getDim() == rhs.getDim()) {
            vector<uint> fValuesL;
            vector<uint> fValuesR;

            for (int i = 0; i <= lhs.getDim(); i++) {
                fValuesL.push_back(filtration[lhs.getConstVertices()[i]]);
                fValuesR.push_back(filtration[rhs.getConstVertices()[i]]);
            }

            sort(fValuesL.begin(), fValuesL.end(), std::greater<uint>());
            sort(fValuesR.begin(), fValuesR.end(), std::greater<uint>());

            for (uint i = 0; i < fValuesL.size(); i++) {
                if (fValuesL[i] == fValuesR[i]) continue;

                return fValuesL[i] < fValuesR[i];
            }
        }

        return lhs.getDim() < rhs.getDim();
    }

    // Compare the filtration index of the two simplexes (vector-valued
    // filtration)
    bool cmpInjectiveFiltr(const implicitS &lhs, const implicitS &rhs) {
        if (lhs.getDim() == rhs.getDim()) {
            vector<uint> fValuesL = simplexFiltration(lhs);
            vector<uint> fValuesR = simplexFiltration(rhs);

            int greater = 0;
            int smaller = 0;

            for (uint i = 0; i < fValuesL.size(); i++) {
                if (fValuesL[i] > fValuesR[i]) greater++;
                if (fValuesL[i] < fValuesR[i]) smaller++;
            }

            if (greater == 0)
                return true;
            else if (greater == smaller) {
                for (uint i = 0; i < fValuesL.size(); i++) {
                    return fValuesL[i] < fValuesR[i];
                }
            } else {
                return false;
            }
        }

        return lhs.getDim() < rhs.getDim();
    }

    bool sortVerticesFiltration(const int &v1, const int &v2) { return filtration[v1] > filtration[v2]; }

public:
    void saveScalarFieldVTK(const char *filename);

    void saveScalarFieldVTK_OG(const char *filename);

    inline void outputFullData(const string &outfile) {
        cout << "\noutput full data\n";
        //    cout << "The Forman gradient has critical cells: " << endl;
        //    for (auto lvl : criticalS) {
        //      cout << "Dim: " << lvl.first << "  #: " << lvl.second.size() << endl;
        //    }
        //    cout << endl;

        auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this, boost::placeholders::_1, boost::placeholders::_2);

        // ciccio's code
        //    SSet cp = SSet(foo);
        //    for (auto lvl : criticalS) {
        //      cp.insert(lvl.second.begin(), lvl.second.end());
        //    }

        // outCriticalPoints(cp);
        // outCriticalPoints_01(cp);
        // saveScalarFieldVTK(outfile.c_str());  // "scalarField.vtk"
        saveScalarFieldVTK_OG((outfile.substr(0, outfile.size() - 4) + "_og.vtk").c_str());
    }

    bool filtrComparer(const pair<float, uint> &v1, const pair<float, uint> &v2) const {
        if (v1.first == v2.first) {
            return v2 < v1;
        }

        return v1.first < v2.first;
    }
};

#endif  // FORMANGRADIENT_H
