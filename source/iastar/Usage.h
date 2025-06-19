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
