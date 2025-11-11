#ifndef USAGE_H
#define USAGE_H

#include <iostream>

#include "fstream"
#include "sstream"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/resource.h"
#include "sys/types.h"
#include "unistd.h"
using namespace std;

// Get current memory usage in KB (Linux-specific, reads from /proc/self/status)
// Returns VmRSS (Resident Set Size - physical memory currently used)
inline long getMemoryUsageKB() {
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

class MemoryUsage {
 private:
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;

 public:
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

#endif  // USAGE_H
