#ifndef USAGE2_H
#define USAGE2_H

#include <iostream>

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/resource.h"
#include "sys/types.h"

using namespace std;

// https://github.com/topology-tool-kit/ttk/blob/591c27277e60f7a6da4c7cceacaa13a20a953108/core/base/common/Os.cpp
// https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
float xx_getMemoryInstantUsage(const bool &output = false) {
    // horrible hack since getrusage() doesn't seem to work well under linux

    // get Virtual Memory (not yet), Physical Memory?
    std::stringstream procFileName;
    // procFileName << "/proc/" << getpid() << "/statm";

    procFileName << "/proc/" << getpid() << "/status";
    std::ifstream procFile(procFileName.str().data(), std::ios::in);
    std::string line;
    if (procFile) {
        float memoryUsage;  // kB
        // procFile >> memoryUsage;
        while (std::getline(procFile, line)) {
            std::istringstream iss(line);
            string str;
            iss >> str;
            if (str.find("VmRSS") != std::string::npos) {
                // VmRSS: 26527016 kB
                iss >> memoryUsage;
            }
        }
        procFile.close();

        if (output) {
            std::cout << "pid: " << getpid() << ", " << memoryUsage / 1024.0 << " MB" << std::endl;
        }
        return memoryUsage / 1024.0;
    }
    return 0;
}

class xx_MemoryUsage {
private:
    int who = RUSAGE_SELF;
    struct rusage usage;
    int ret;

public:
    // old code from Ciccio. Seems not work.
    inline float getValue_in_KB(bool output) {  // Note: this value is in KB!

        ret = getrusage(who, &usage);
        if (output) cout << "Memory Usage: " << usage.ru_maxrss / (1024.0) << " KB" << endl;

        return usage.ru_maxrss / (1024.0);
    }

    inline float getValue_in_MB(bool output) {  // Note: this value is in MB!

        ret = getrusage(who, &usage);
        if (output) cout << "Memory Usage: " << usage.ru_maxrss / (1024.0 * 1024.0) << " MB" << endl;

        return usage.ru_maxrss / (1024.0 * 1024.0);
    }

    inline float getValue_in_GB(bool output) {  // Note: this value is in GB!

        ret = getrusage(who, &usage);
        if (output) cout << "Memory Usage: " << usage.ru_maxrss / (1024.0 * 1024.0 * 1024.0) << " GB" << endl;

        return usage.ru_maxrss / (1024.0 * 1024.0 * 1024.0);
    }
};

#endif  // USAGE2_H
