#include <numeric>

#include "formangradient.h"

// original
// Compute DescendingCells (desCells) for the given cell (cell).
// given cell should be a critical simplex. >=1
// DO NOT MODFIY the code.
void FormanGradient::computeDescendingCell(bool output, implicitS const &cell, SSet &desCells) {
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this, boost::placeholders::_1, boost::placeholders::_2);
    SSet cellsPath = SSet(foo);
    if (output) {
        desCells = SSet(foo);
        cellsPath.insert(cell);
    }
    stack<implicitS> qu;
    qu.push(cell);

    while (!qu.empty()) {
        implicitS top = qu.top();
        qu.pop();

        vector<implicitS> *bd = sc.boundaryk(top, top.getDim() - 1);

        for (auto s: *bd) {
            implicitS next;
            if (getPair(s, next)) {
                assert(isPaired(next) && isPaired(s));
                if (next.getDim() == cell.getDim() && !sc.theSame(next, top)) {
                    qu.push(next);
                    if (output) {
                        cellsPath.insert(next);
                    }
                }
            } else {
                // meet critical cells. do nothing
                //        if (output) {
                //          desCells.insert(cellsPath.begin(), cellsPath.end());
                //          cellsPath.clear();
                //        }

                // original
                //        if (cell.getDim() == 1 /* && validConnection(cell,s)*/) {
                //          desCells.insert(cellsPath.begin(), cellsPath.end());
                //          cellsPath.clear();
                //        }
            }
        }
        delete bd;

        // finish visiting all boundaries. add to desCells
        if (output) {
            desCells.insert(cellsPath.begin(), cellsPath.end());
            cellsPath.clear();
        }
    }
}

/// original
// Compute AscendingCell (ascCells) for the given cell (cell).
// given cell should be a critical simplex. not the max simplex in the data

// warning! used in the TOPOSEGMENT. do not modify!
// follow Ciccio book
/*
 * 2D case
 * a descending 2-cell corresponds to a maximum: a collection of triangles
 * an ascending 2-cell on minimum: a collection of vertices
 * a descending 1-cell on saddle: a collection of edges
 * an ascending 1-cell on saddle: a collection of edges
 */
void FormanGradient::computeAscendingCell(bool output, implicitS const &cell, SSet &ascCells) {
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this, boost::placeholders::_1, boost::placeholders::_2);
    SSet cellsPath = SSet(foo);
    if (output) {
        ascCells = SSet(foo);
        cellsPath.insert(cell);
    }
    queue<implicitS> qu;
    qu.push(cell);
    while (!qu.empty()) {
        //        cout << "Navigo" << endl;
        implicitS top = qu.front();
        qu.pop();
        vector<implicitS> *bd = sc.coboundaryk(top, top.getDim() + 1);

        for (auto s: *bd) {
            // *bd
            implicitS next;
            if (getPair(s, next)) {
                if (next.getDim() == cell.getDim() && !sc.theSame(next, top)) {
                    qu.push(next);
                    if (output) {
                        cellsPath.insert(next);
                    }
                } else {
                }
            } else {
                //        cout << "is critical edge\n";
                //        if (cell.getDim() == 0 /*&& validConnection(cell,s)*/) {
                //          ascCells.insert(cellsPath.begin(), cellsPath.end());
                //          cellsPath.clear();
                //        }
            }
        }
        delete bd;
        // finish visiting all co-boundary
        // add to ascCells
        if (output) {
            ascCells.insert(cellsPath.begin(), cellsPath.end());
            cellsPath.clear();
        }
    }
}


void FormanGradient::computeAscendingCell_xx(bool output, implicitS const &cell, SSet &ascCells) {
    // diff from non-xx version.
    /* non-xx version: input vertex critical cell, output: vertex
     * xx version: input vertex critical cell, output: saddle.  input saddle critical cells, output: triangles
     */
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this, boost::placeholders::_1, boost::placeholders::_2);
    SSet cellsPath = SSet(foo);
    if (output) {
        ascCells = SSet(foo);
        cellsPath.insert(cell);
    }
    queue<implicitS> qu;
    qu.push(cell);
    while (!qu.empty()) {
        //        cout << "Navigo" << endl;
        implicitS top = qu.front();
        qu.pop();
        vector<implicitS> *bd = sc.coboundaryk(top, top.getDim() + 1);

        for (auto s: *bd) {
            // *bd
            implicitS next;
            if (getPair(s, next)) {
                if (next.getDim() == cell.getDim() && !sc.theSame(next, top)) {
                    qu.push(next);
                    if (output) {
                        cellsPath.insert(s);
                    }
                } else {
                }
            } else {
                //        cout << "is critical edge\n";
                //if (cell.getDim() == 0 /*&& validConnection(cell,s)*/) {
                if (output) {
                    ascCells.insert(cellsPath.begin(), cellsPath.end());
                    cellsPath.clear();
                    //}
                }
            }
        }
        delete bd;
        // finish visiting all co-boundary
        // add to ascCells
        //        if (output) {
        //            ascCells.insert(cellsPath.begin(), cellsPath.end());
        //            cellsPath.clear();
        //        }
    }
}

/// visualization:
// Ciccio
void FormanGradient::visMorse() {
    list<SSet> descending1manifold;
    list<SSet> ascending1manifold;
    for (auto criticalLVL: criticalS) {
        for (implicitS c: criticalLVL.second) {
            SSet ascending, descending;

            if (c.getDim() == 1) {
                computeDescendingCell(true, c, descending);
                descending1manifold.push_back(descending);

                //computeAscendingCell_xx(true, c, ascending);
                computeAscendingCell(true, c, ascending);
                ascending1manifold.push_back(ascending);
            }

            // xx's test
            //      SSet descending2;
            //      if (c.getDim() == 3) {
            //        computeDescendingCell(true, c, descending2);
            //        descending2manifold.push_back(descending2);
            //      }
        }
    }

    out1cells(descending1manifold, true);
    out1cells(ascending1manifold, false);
}
