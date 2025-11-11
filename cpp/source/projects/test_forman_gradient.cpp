//
// Simple test program for Forman gradient computation
// Tests only the essential steps: reading, building data structure, and computing Forman gradient
// Outputs timing and file statistics for performance analysis
//

#include <iostream>
#include <string>
#include "../forman/formangradient.h"
#include "../iastar/Timer.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input_file> [funID]" << endl;
        cout << "  input_file: .off or .ply file" << endl;
        cout << "  funID: scalar function ID (default: 3)" << endl;
        return 1;
    }
    
    string infile = argv[1];
    int funID = 3; // default
    if (argc >= 3) {
        funID = atoi(argv[2]);
    }
    
    cout << "========================================" << endl;
    cout << "Forman Gradient Test Program" << endl;
    cout << "========================================" << endl;
    cout << "Input file: " << infile << endl;
    cout << "Function ID: " << funID << endl;
    cout << endl;
    
    IO_Timer total_timer;
    total_timer.start();
    
    // Create FormanGradient - this will internally time reading and gradient computation
    FormanGradient fg(infile, funID);
    
    // Get file statistics
    cout << "\n=== File Statistics ===" << endl;
    cout << "Vertices: " << fg.sc.getVerticesNum() << endl;
    
    size_t total_top_simplexes = 0;
    auto top_set = fg.sc.getTopSimplexesSet();
    cout << "Top Simplexes by dimension:" << endl;
    for (int dim : top_set) {
        int count = fg.sc.getTopSimplexesNum(dim);
        total_top_simplexes += count;
        cout << "  Dim " << dim << ": " << count << endl;
    }
    cout << "Total top simplexes: " << total_top_simplexes << endl;
    cout << "Complex dimension: " << fg.sc.getComplexDim() << endl;
    
    total_timer.stop();
    cout << "\n=== Total Execution Time ===" << endl;
    cout << "Total time: " << total_timer.getElapsedTime() << " s" << endl;
    
    cout << "\nTest completed successfully!" << endl;
    
    return 0;
}

