//
// Simple test program for Forman gradient computation
// Tests only the essential steps: reading, building data structure, and computing Forman gradient
// Outputs timing and file statistics for performance analysis
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "../forman/formangradient.h"
#include "../iastar/Timer.h"

using namespace std;

// Function to get current memory usage in KB (Linux-specific)
long getMemoryUsageKB() {
    long memory_kb = 0;
    ifstream status_file("/proc/self/status");
    if (status_file.is_open()) {
        string line;
        while (getline(status_file, line)) {
            if (line.find("VmRSS:") == 0) {  // Resident Set Size (physical memory)
                istringstream iss(line);
                string key, value, unit;
                iss >> key >> value >> unit;
                memory_kb = stol(value);
                break;
            }
        }
        status_file.close();
    }
    return memory_kb;
}

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
    
    // Measure initial memory usage
    long memory_before = getMemoryUsageKB();
    
    // Create FormanGradient - this will internally time reading and gradient computation
    FormanGradient fg(infile, funID);
    
    // Measure memory usage after building data structures
    long memory_after = getMemoryUsageKB();
    long memory_used = memory_after - memory_before;
    
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
    
    cout << "\n=== Memory Usage ===" << endl;
    cout << "Memory before: " << memory_before << " KB" << endl;
    cout << "Memory after: " << memory_after << " KB" << endl;
    cout << "Memory used: " << memory_used << " KB" << endl;
    cout << "Memory used (MB): " << (memory_used / 1024.0) << " MB" << endl;
    
    // Output in parseable format for performance report
    cout << "Memory usage (KB): " << memory_after << endl;
    cout << "Memory usage (MB): " << (memory_after / 1024.0) << endl;
    
    cout << "\nTest completed successfully!" << endl;
    
    return 0;
}

