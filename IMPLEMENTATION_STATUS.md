# Performance Testing Implementation Status

## ✅ Completed Components

### 1. Comparison Tool (`cmp_ply.py`)
- **Location**: `python/tests/cmp_ply.py`
- **Status**: ✅ Complete
- **Features**:
  - Compares two `.ply` files with segmentation labels
  - Computes IoU (Intersection over Union) per label and overall
  - Reports detailed comparison results
  - Exit code: 0 if identical, 1 if different

**Usage**:
```bash
. ~/xx_pyvenvs/treemapping_project/bin/activate.fish
python python/tests/cmp_ply.py file1.ply file2.ply
```

### 2. Enhanced Timing Instrumentation
- **Status**: ✅ Complete
- **Files Modified**:
  - `cpp/source/iastar/io_functions.cpp`:
    - `readPLY()`: Enhanced timing for file I/O, vertex parsing, cell parsing, build time, total time
    - `readOFF()`: Enhanced timing for file I/O, vertex parsing, cell parsing, build time, total time
    - Added `USE_PARALLEL_BUILD` compile-time flag support
  - `cpp/source/forman/formangradient.cpp`:
    - Enhanced timing for gradient encoding and filtration computation
    - Separate timers for each major phase

**Timing Output Format**:
```
=== Reading Input ===
   [off/ply] read file I/O: X.XXX s
   [off/ply] read vertices: X.XXX s
   [off/ply] read cells: X.XXX s
   [off/ply] build IA* (sequential/parallel): X.XXX s
   [off/ply] total reading time: X.XXX s

=== Computing Forman Gradient ===
   Gradient encoding time: X.XXX s
   Filtration computation time: X.XXX s
   Total Forman gradient time: X.XXX s
```

### 3. Test Program (`test_forman_gradient.cpp`)
- **Location**: `cpp/source/projects/test_forman_gradient.cpp`
- **Status**: ✅ Complete
- **Features**:
  - Lightweight program for testing Forman gradient computation only
  - Outputs file statistics (vertex count, top simplex count)
  - Outputs detailed timing information
  - No full segmentation workflow

**Usage**:
```bash
./test_forman_gradient <input_file> [funID]
```

### 4. CMake Configuration
- **Status**: ✅ Complete
- **Files Modified**:
  - `cpp/CMakeLists.txt`: Added `test_forman_gradient` executable target
- **Compile-time Flags**:
  - `USE_PARALLEL_BUILD`: Controls sequential vs parallel data structure building
    - `ON`: Uses `buildDataStructure_parallel()`
    - `OFF` (default): Uses `buildDataStructure()`

### 5. Compilation Script
- **Location**: `compile_variants.sh`
- **Status**: ✅ Complete
- **Features**:
  - Automates compilation of 4 executable variants
  - Handles git branch switching for old/new code
  - Creates executables:
    - `xx_tts_seq_old`: Old code, sequential
    - `xx_tts_pa_old`: Old code, parallel
    - `xx_tts_seq_new`: New code, sequential
    - `xx_tts_pa_new`: New code, parallel

**Usage**:
```bash
./compile_variants.sh [build_dir]
```

---

## ⚠️ Pending Tasks

### 1. Git Branch Setup
- **Status**: ⚠️ Needs manual action
- **Action Required**:
  ```bash
  # Find commit before adjRelations removal
  git log --oneline --grep="adjRelations" -10
  
  # Create old branch (replace <commit-hash> with actual hash)
  git checkout -b old_adjRelations <commit-hash>
  
  # Return to main/new branch
  git checkout main  # or your current branch
  ```

### 2. Compilation
- **Status**: ⚠️ Ready to compile
- **Steps**:
  1. Set up git branch (see above)
  2. Run compilation script:
     ```bash
     ./compile_variants.sh cpp/cmake-build-debug
     ```
  3. Verify executables are created:
     ```bash
     ls -lh cpp/cmake-build-debug/xx_tts_*
     ```

### 3. Testing
- **Status**: ⚠️ Ready to test
- **Test Files**:
  - `close_stems_3_a0.01.ply` - Complete workflow (segmentation)
  - `close_stems_3_a0.010.off` - Forman gradient only
  - `aoi_thin_low_pts_a0.010.off` - Forman gradient only
  - `tree_228_veg_xyz_as_0.01.off` - Forman gradient only

- **Test Steps**:
  1. Run complete workflow on `.ply` file with all 4 executables
  2. Compare segmentation results using `cmp_ply.py`
  3. Run Forman gradient tests on all `.off` files
  4. Collect timing data
  5. Analyze scalability (file size vs performance)

---

## 📋 Next Steps

1. **Create git branch for old code**:
   ```bash
   git log --oneline | head -20  # Find commit before changes
   git checkout -b old_adjRelations <commit-hash>
   git checkout main  # Return to current branch
   ```

2. **Compile all variants**:
   ```bash
   ./compile_variants.sh
   ```

3. **Run tests**:
   - For `.ply` file: Run complete workflow, compare results
   - For `.off` files: Run Forman gradient tests, collect timing

4. **Collect and analyze data**:
   - Parse timing output
   - Create scalability tables
   - Compare old vs new performance
   - Document findings

---

## 🔧 Technical Notes

### Compile-Time Flags

The code uses `USE_PARALLEL_BUILD` preprocessor flag to control sequential vs parallel builds:

```cpp
#ifdef USE_PARALLEL_BUILD
    buildDataStructure_parallel();
#else
    buildDataStructure();
#endif
```

To set this flag in CMake:
```cmake
# For parallel build
cmake .. -DUSE_PARALLEL_BUILD=ON

# For sequential build (default)
cmake .. -DUSE_PARALLEL_BUILD=OFF
```

### File Statistics

File statistics (vertex count, top simplex count) are automatically printed by:
- `readOFF()` and `readPLY()` functions
- `test_forman_gradient` program

These can be parsed from output logs for scalability analysis.

---

## 📊 Expected Output Format

### From `test_forman_gradient`:
```
========================================
Forman Gradient Test Program
========================================
Input file: close_stems_3_a0.010.off
Function ID: 3

=== Reading Input ===
   off read file I/O: X.XXX s
   off read vertices: X.XXX s
   off read cells: X.XXX s
   off build IA* (sequential): X.XXX s
   off total reading time: X.XXX s

=== Computing Forman Gradient ===
   Gradient encoding time: X.XXX s
   Filtration computation time: X.XXX s
   Total Forman gradient time: X.XXX s

=== File Statistics ===
Vertices: 12345
Top Simplexes by dimension:
  Dim 2: 67890
Total top simplexes: 67890
Complex dimension: 2

=== Total Execution Time ===
Total time: X.XXX s
```

### From `cmp_ply.py`:
```
============================================================
PLY Segmentation Comparison Tool
============================================================

Reading file1.ply...
Reading file2.ply...
Comparing 12345 points...

Per-label IoU:
------------------------------------------------------------
Label      IoU            Count (File1)   Count (File2)   
------------------------------------------------------------
0          1.000000       1000            1000            ✓
1          1.000000       2000            2000            ✓
------------------------------------------------------------

Overall IoU: 1.000000
✓ Files are IDENTICAL
```

---

**Last Updated**: 2025-01-XX
**Status**: Implementation Complete, Ready for Testing

