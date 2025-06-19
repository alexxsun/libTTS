#include <iostream>
#include <iostream>
#include "projects/TopoAnalysis.h"

int main(int argc, char **argv) {
    std::string infile = argv[1];
    TopoAnalysis ta(infile, 3, true);
    return 0;
}
