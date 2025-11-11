# Performance Testing and Validation Plan

## Overview

This plan outlines the comprehensive testing strategy to validate the performance improvements from removing `adjRelations` and using `completeCoboundaryTop`. The plan includes:

1. **Segmentation Comparison Tool** (`cmp_ply.py`) - Validate correctness
2. **Enhanced Timing Instrumentation** - Track key performance metrics
3. **Simple Forman Gradient Test Program** - Isolated performance testing
4. **Four Executable Variants** - Compare old vs new, sequential vs parallel
5. **Comprehensive Test Suite** - Multiple test files with different workflows

---

## 1. Segmentation Comparison Tool (`cmp_ply.py`)

### Purpose
Compare two `.ply` segmentation files to ensure the new code produces identical results to the old code.

### Implementation Details

**Location**: `/home/alex/Projects/libTTS_public/python/tests/cmp_ply.py`

**Functionality**:
- Read two `.ply` files with `x`, `y`, `z`, `label` properties
- Extract vertex data: `vertex_data = np.vstack([plydata['vertex']['x'], plydata['vertex']['y'], plydata['vertex']['z'], plydata['vertex']['label']]).T`
- Compare segmentations using **Intersection over Union (IoU)** per label
- Overall IoU should be **1.0** if files are identical

**Key Features**:
- Use `plyfile` library (available in virtual environment)
- Handle label mismatches gracefully
- Report per-label IoU and overall IoU
- Exit code: 0 if identical, 1 if different

**Usage**:
```bash
. ~/xx_pyvenvs/treemapping_project/bin/activate.fish
python python/tests/cmp_ply.py file1.ply file2.ply
```

**Output Format**:
```
Comparing file1.ply and file2.ply...
Label 0: IoU = 1.000 (100% match)
Label 1: IoU = 1.000 (100% match)
...
Overall IoU: 1.000
Files are IDENTICAL ✓
```

---

## 2. Enhanced Timing Instrumentation

### Purpose
Track and display timing for critical operations to measure performance improvements.

### Key Metrics to Track

1. **Total Reading Time** (`readOFF`/`readPLY`)
   - File I/O time
   - Vertex parsing time
   - Top simplex parsing time
   - **Data structure building time** (IA* build)

2. **Data Structure Build Time** (IA*)
   - Sequential build time (`buildDataStructure()`)
   - Parallel build time (`buildDataStructure_parallel()`)
   - Breakdown by phase (if possible)

3. **Tops Computation Time**
   - Time to compute topological features
   - (If applicable in the workflow)

4. **Forman Gradient Computation Time**
   - Gradient encoding creation
   - Filtration computation
   - Component-based filtration

5. **Memory Usage** (optional but recommended)
   - Peak memory during data structure building
   - Peak memory during Forman gradient computation

### Implementation Changes

#### Update `readOFF(const char *file, const int &funID)`

**Current State**:
- Has some timing but inconsistent with `readPLY`
- Uses `buildDataStructure()` or `buildDataStructure_parallel()` but not consistently timed

**Required Changes**:
```cpp
void SimplicialComplex::readOFF(const char *file, const int &funID) {
    IO_Timer total_timer, io_timer, build_timer;
    
    total_timer.start();
    
    // File I/O and parsing
    io_timer.start();
    // ... existing file reading code ...
    io_timer.stop();
    cout << "   off read file I/O: " << io_timer.getElapsedTime() << " s" << endl;
    
    // Data structure building
    build_timer.start();
    buildDataStructure();  // or buildDataStructure_parallel() based on flag
    build_timer.stop();
    cout << "   off build IA* (seq): " << build_timer.getElapsedTime() << " s" << endl;
    // OR
    cout << "   off build IA* (parallel): " << build_timer.getElapsedTime() << " s" << endl;
    
    total_timer.stop();
    cout << "   off total reading time: " << total_timer.getElapsedTime() << " s" << endl;
}
```

#### Update `readPLY(const char *file)`

**Current State**:
- Has timing for vertices and cells
- Has timing for build but inconsistent format

**Required Changes**:
- Make timing format consistent with `readOFF`
- Add total reading time
- Clearly separate sequential vs parallel build timing

#### Update Forman Gradient Computation

**Location**: `cpp/source/forman/formangradient.cpp`

**Current State**:
- Has some timing but not comprehensive

**Required Changes**:
```cpp
// In FormanGradient constructor or compute method
IO_Timer forman_timer;

forman_timer.start();
// Gradient encoding creation
gradient = GradientEncoding(sc);
forman_timer.stop();
cout << "Forman gradient encoding time: " << forman_timer.getElapsedTime() << " s" << endl;

forman_timer.start();
// Filtration computation
// ... filtration code ...
forman_timer.stop();
cout << "Filtration computation time: " << forman_timer.getElapsedTime() << " s" << endl;
```

### Timing Output Format

**Standardized Format**:
```
=== Reading Input ===
   [off/ply] read file I/O: X.XXX s
   [off/ply] read vertices: X.XXX s
   [off/ply] read cells: X.XXX s
   [off/ply] build IA* (seq/parallel): X.XXX s
   [off/ply] total reading time: X.XXX s

=== Forman Gradient Computation ===
   Gradient encoding time: X.XXX s
   Filtration computation time: X.XXX s
   Total Forman gradient time: X.XXX s
```

---

## 3. Simple Forman Gradient Test Program

### Purpose
Create a lightweight test program that only runs the essential steps needed to test Forman gradient computation, without the full segmentation workflow.

### Implementation

**Location**: `cpp/source/projects/test_forman_gradient.cpp`

**Functionality**:
1. Read input file (`.off` or `.ply`)
2. Build data structure (IA*)
3. Compute Forman gradient
4. Output timing and memory usage
5. Exit (no segmentation)

**Key Features**:
- Minimal dependencies
- Clear timing output
- Optional memory profiling
- Can be compiled with different flags (old/new, seq/parallel)

**Code Structure**:
```cpp
#include "formangradient.h"
#include <iostream>
#include <chrono>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }
    
    string infile = argv[1];
    int funID = 3; // default
    
    cout << "=== Forman Gradient Test ===" << endl;
    cout << "Input file: " << infile << endl;
    
    IO_Timer total_timer;
    total_timer.start();
    
    // This will internally time reading and building
    FormanGradient fg(infile, funID);
    
    // Additional timing for gradient computation if needed
    // (FormanGradient constructor may already handle this)
    
    total_timer.stop();
    cout << "\n=== Total Time ===" << endl;
    cout << "Total execution time: " << total_timer.getElapsedTime() << " s" << endl;
    
    return 0;
}
```

**CMake Integration**:
- Add to `cpp/CMakeLists.txt` or `cpp/source/CMakeLists.txt`
- Create separate targets for old/new versions

---

## 4. Four Executable Variants

### Purpose
Compile four versions to compare:
- **Old code** (with `adjRelations`) vs **New code** (with `completeCoboundaryTop`)
- **Sequential** (`buildDataStructure()`) vs **Parallel** (`buildDataStructure_parallel()`)

### Executables to Create

1. **`xx_tts_seq_old`**: Old code, sequential build
2. **`xx_tts_pa_old`**: Old code, parallel build
3. **`xx_tts_seq_new`**: New code, sequential build
4. **`xx_tts_pa_new`**: New code, parallel build

### Implementation Strategy

#### Option A: Git Branches (Recommended)
1. Create a branch `old_adjRelations` from current state (before removal)
2. Keep `main` branch with new code
3. Compile from both branches with different flags

#### Option B: Compile-Time Flags
1. Add `#ifdef USE_ADJ_RELATIONS` guards
2. Compile with/without the flag
3. More complex but single codebase

#### Option C: Separate Directories
1. Copy code to `cpp_old/` and `cpp_new/`
2. Compile separately
3. Simple but code duplication

**Recommended: Option A (Git Branches)**

### Compilation Steps

```bash
# 1. Save current state as "new" branch
git checkout -b new_completeCoboundaryTop
git add -A
git commit -m "New code with completeCoboundaryTop"

# 2. Create "old" branch from commit before adjRelations removal
git checkout -b old_adjRelations <commit-before-removal>

# 3. Compile old versions
git checkout old_adjRelations
cd cpp/cmake-build-debug
# Modify readOFF/readPLY to use buildDataStructure() for seq
cmake .. -DCMAKE_BUILD_TYPE=Release
make xx_tts_seq_old
# Modify readOFF/readPLY to use buildDataStructure_parallel() for parallel
make xx_tts_pa_old

# 4. Compile new versions
git checkout new_completeCoboundaryTop
cd cpp/cmake-build-debug
# Modify readOFF/readPLY to use buildDataStructure() for seq
cmake .. -DCMAKE_BUILD_TYPE=Release
make xx_tts_seq_new
# Modify readOFF/readPLY to use buildDataStructure_parallel() for parallel
make xx_tts_pa_new
```

### Code Modifications Needed

**For Sequential Versions**:
- In `readOFF()` and `readPLY()`, ensure `buildDataStructure()` is called (not `buildDataStructure_parallel()`)

**For Parallel Versions**:
- In `readOFF()` and `readPLY()`, ensure `buildDataStructure_parallel()` is called

**Current State Check**:
- `readOFF(const char *file, const int &funID)`: Uses `buildDataStructure_parallel()` (line 375)
- `readPLY(const char *file)`: Uses `buildDataStructure()` (line 204)

**Action Required**:
- Add a compile-time flag or runtime parameter to choose sequential vs parallel
- OR create wrapper functions
- OR modify directly for each executable variant

---

## 5. Test Files and Workflows

### Test Files Location
`/home/alex/Projects/libTTS_public/cpp/cmake-build-debug/`

### Test Files

1. **`close_stems_3_a0.01.ply`**
   - **Workflow**: Complete workflow (segmentation)
   - **Purpose**: Validate correctness - compare segmentation results
   - **Comparison**: Run all 4 executables, compare output `.ply` files using `cmp_ply.py`
   - **Expected**: All 4 should produce identical segmentation (IoU = 1.0)

2. **`close_stems_3_a0.010.off`**
   - **Workflow**: Forman gradient only (no segmentation)
   - **Purpose**: Measure performance for data structure building and Forman gradient
   - **Comparison**: Compare timing across 4 executables

3. **`aoi_thin_low_pts_a0.010.off`**
   - **Workflow**: Forman gradient only
   - **Purpose**: Measure performance on different dataset

4. **`tree_228_veg_xyz_as_0.01.off`**
   - **Workflow**: Forman gradient only
   - **Purpose**: Measure performance on different dataset

### Test Execution Plan

#### For `.ply` file (Complete Workflow)

```bash
cd /home/alex/Projects/libTTS_public/cpp/cmake-build-debug

# Run all 4 versions
./xx_tts_seq_old <args_for_complete_workflow> close_stems_3_a0.01.ply
./xx_tts_pa_old <args_for_complete_workflow> close_stems_3_a0.01.ply
./xx_tts_seq_new <args_for_complete_workflow> close_stems_3_a0.01.ply
./xx_tts_pa_new <args_for_complete_workflow> close_stems_3_a0.01.ply

# Compare results
. ~/xx_pyvenvs/treemapping_project/bin/activate.fish
python ../../python/tests/cmp_ply.py output_seq_old.ply output_seq_new.ply
python ../../python/tests/cmp_ply.py output_pa_old.ply output_pa_new.ply
python ../../python/tests/cmp_ply.py output_seq_old.ply output_pa_old.ply  # Should be same
python ../../python/tests/cmp_ply.py output_seq_new.ply output_pa_new.ply  # Should be same
```

#### For `.off` files (Forman Gradient Only)

```bash
cd /home/alex/Projects/libTTS_public/cpp/cmake-build-debug

# Option 1: Use test_forman_gradient program
./test_forman_gradient_seq_old close_stems_3_a0.010.off > timing_seq_old.log
./test_forman_gradient_pa_old close_stems_3_a0.010.off > timing_pa_old.log
./test_forman_gradient_seq_new close_stems_3_a0.010.off > timing_seq_new.log
./test_forman_gradient_pa_new close_stems_3_a0.010.off > timing_pa_new.log

# Option 2: Use xx_tts with minimal workflow (if test program not ready)
# (Need to check what minimal args are needed)
```

### Expected Results

#### Correctness (`.ply` file)
- All 4 executables should produce **identical** segmentation results
- IoU between any two outputs should be **1.0**

#### Performance (`.off` files)
- **New code should be faster** than old code (especially for data structure building)
- **Parallel should be faster** than sequential (for large datasets)
- **Memory usage** may be higher for new code (storing complete mapping)

### Performance Metrics to Collect

For each test file and each executable:

1. **Reading Time Breakdown**:
   - File I/O time
   - Vertex parsing time
   - Top simplex parsing time
   - Data structure build time (IA*)

2. **Forman Gradient Time**:
   - Gradient encoding time
   - Filtration computation time
   - Total Forman gradient time

3. **Total Execution Time**

4. **Memory Usage** (if available):
   - Peak memory during data structure building
   - Peak memory during Forman gradient computation

5. **File Statistics** (for scalability analysis):
   - Number of vertices (points)
   - Number of top simplexes (by dimension)
   - File size on disk (optional)

---

## 6. Scalability Analysis: File Size vs Performance

### Purpose
Understand how performance scales with dataset size by analyzing the relationship between point count and execution time/memory usage across the three test files.

### Data Collection

For each test file, collect:

1. **File Characteristics**:
   - Number of vertices (`getVerticesNum()`)
   - Number of top simplexes (total and by dimension)
   - File size on disk (optional, for reference)

2. **Performance Metrics** (from each executable):
   - Data structure build time (IA*)
   - Forman gradient computation time
   - Total execution time
   - Peak memory usage (if available)

### Implementation

#### Option 1: Automatic Collection Script

Create a Python script `collect_performance_data.py` that:
- Parses timing output from test runs
- Extracts file statistics from output logs
- Generates summary tables and plots

#### Option 2: Enhanced Test Program Output

Modify `test_forman_gradient.cpp` to output structured data:

```cpp
// Output format (CSV or structured text):
// File: close_stems_3_a0.010.off
// Vertices: 12345
// Top Simplexes: 67890
// Build IA* Time: 1.234 s
// Forman Gradient Time: 0.567 s
// Total Time: 1.801 s
// Peak Memory: 123.45 MB
```

#### Option 3: Manual Collection Template

Create a spreadsheet template with columns:
- File name
- # Vertices
- # Top Simplexes
- Build Time (seq_old, pa_old, seq_new, pa_new)
- Forman Gradient Time (seq_old, pa_old, seq_new, pa_new)
- Total Time (seq_old, pa_old, seq_new, pa_new)
- Memory Usage (seq_old, pa_old, seq_new, pa_new)

### Analysis Output

#### Summary Table Format

| File | # Vertices | # Top Simplexes | Build IA* (seq_new) | Build IA* (pa_new) | Forman (seq_new) | Forman (pa_new) | Total (seq_new) | Total (pa_new) | Memory (new) |
|------|------------|-----------------|---------------------|---------------------|------------------|-----------------|-----------------|----------------|--------------|
| close_stems_3_a0.010.off | X | Y | T1 | T2 | T3 | T4 | T5 | T6 | M1 |
| aoi_thin_low_pts_a0.010.off | X | Y | T1 | T2 | T3 | T4 | T5 | T6 | M1 |
| tree_228_veg_xyz_as_0.01.off | X | Y | T1 | T2 | T3 | T4 | T5 | T6 | M1 |

#### Performance per Vertex Metrics

Calculate and display:
- **Time per vertex** (build time / # vertices)
- **Time per top simplex** (build time / # top simplexes)
- **Memory per vertex** (peak memory / # vertices)

| File | # Vertices | Build Time/Vertex (ms) | Forman Time/Vertex (ms) | Memory/Vertex (KB) |
|------|------------|------------------------|-------------------------|---------------------|
| close_stems_3_a0.010.off | X | Y1 | Y2 | Y3 |
| aoi_thin_low_pts_a0.010.off | X | Y1 | Y2 | Y3 |
| tree_228_veg_xyz_as_0.01.off | X | Y1 | Y2 | Y3 |

#### Scalability Insights

Look for:
1. **Linear scaling**: Time increases proportionally with vertex count
2. **Sub-linear scaling**: Time increases slower than vertex count (good!)
3. **Super-linear scaling**: Time increases faster than vertex count (may indicate bottlenecks)
4. **Memory scaling**: How memory usage grows with dataset size

#### Comparison Matrix: Old vs New

| Metric | seq_old | seq_new | Speedup | pa_old | pa_new | Speedup |
|--------|---------|---------|---------|--------|--------|---------|
| **File 1** (small) | | | | | | |
| Build IA* | T1 | T2 | X.Xx | T3 | T4 | X.Xx |
| Forman | T5 | T6 | X.Xx | T7 | T8 | X.Xx |
| Total | T9 | T10 | X.Xx | T11 | T12 | X.Xx |
| **File 2** (medium) | | | | | | |
| Build IA* | T1 | T2 | X.Xx | T3 | T4 | X.Xx |
| Forman | T5 | T6 | X.Xx | T7 | T8 | X.Xx |
| Total | T9 | T10 | X.Xx | T11 | T12 | X.Xx |
| **File 3** (large) | | | | | | |
| Build IA* | T1 | T2 | X.Xx | T3 | T4 | X.Xx |
| Forman | T5 | T6 | X.Xx | T7 | T8 | X.Xx |
| Total | T9 | T10 | X.Xx | T11 | T12 | X.Xx |

### Expected Patterns

1. **Build IA* Time**:
   - Should scale roughly linearly with vertex count
   - New code should show consistent speedup (2-5x) across all file sizes
   - Parallel version should show better speedup on larger files

2. **Forman Gradient Time**:
   - May scale differently depending on topology
   - New code should be similar or slightly faster
   - Less affected by parallelization

3. **Memory Usage**:
   - Should scale roughly linearly with vertex count
   - New code may use 10-30% more memory consistently

4. **Total Time**:
   - Should reflect the sum of components
   - Speedup should be consistent across file sizes

### Visualization (Optional)

Create plots showing:
- **Time vs # Vertices**: Scatter plot with trend lines
- **Memory vs # Vertices**: Scatter plot with trend lines
- **Speedup vs File Size**: Bar chart comparing old vs new

### Implementation Steps

1. **Collect file statistics**:
   - The `readOFF()` and `readPLY()` functions already output vertex and top simplex counts
   - Example output: `"Complex vertices #: 12345"` and `"Complex top simplices #: 67890"`
   - Parse these from output logs or add structured output to test program
   - Alternatively, use `getVerticesNum()` and `getTopSimplexesNum()` in test program

2. **Collect performance data**:
   - Run all 4 executables on all 3 test files
   - Capture timing output to log files

3. **Parse and aggregate**:
   - Use script or manual entry into spreadsheet
   - Calculate per-vertex metrics

4. **Generate summary**:
   - Create comparison tables
   - Calculate speedups
   - Identify scaling patterns

5. **Document findings**:
   - Add to performance report
   - Note any anomalies or unexpected patterns

---

## 7. Implementation Checklist

### Phase 1: Preparation
- [ ] Create `cmp_ply.py` comparison tool
- [ ] Test `cmp_ply.py` on known identical/different files
- [ ] Create git branch for old code (`old_adjRelations`)
- [ ] Verify old code compiles and runs

### Phase 2: Code Updates
- [ ] Update `readOFF()` with enhanced timing
- [ ] Update `readPLY()` with enhanced timing
- [ ] Update Forman gradient timing
- [ ] Create `test_forman_gradient.cpp` program
- [ ] Add CMake targets for test program

### Phase 3: Compilation
- [ ] Compile `xx_tts_seq_old` (old code, sequential)
- [ ] Compile `xx_tts_pa_old` (old code, parallel)
- [ ] Compile `xx_tts_seq_new` (new code, sequential)
- [ ] Compile `xx_tts_pa_new` (new code, parallel)
- [ ] Compile test programs (if using)

### Phase 4: Testing
- [ ] Run complete workflow on `close_stems_3_a0.01.ply` with all 4 executables
- [ ] Compare segmentation results using `cmp_ply.py`
- [ ] Run Forman gradient tests on all `.off` files
- [ ] Collect timing data for all tests
- [ ] Document results

### Phase 5: Analysis
- [ ] Collect file statistics (vertex count, top simplex count) for all test files
- [ ] Compare timing across old vs new
- [ ] Compare timing across sequential vs parallel
- [ ] Analyze memory usage differences
- [ ] Calculate per-vertex performance metrics
- [ ] Analyze scalability patterns (time vs file size)
- [ ] Create scalability summary table
- [ ] Create summary report with scalability insights

---

## 8. Expected Outcomes

### Correctness
- ✅ All 4 executables produce **identical** segmentation results
- ✅ IoU = 1.0 for all comparisons

### Performance Improvements (Expected)
- ✅ **Data structure building**: 2-5x faster (removed expensive cluster expansion)
- ✅ **Forman gradient**: Similar or slightly faster (due to faster `topStar()`)
- ✅ **Total time**: 1.5-3x faster overall
- ⚠️ **Memory**: 10-30% increase (storing complete mapping vs partial)

### Performance Comparison Matrix

| Metric | seq_old | pa_old | seq_new | pa_new |
|--------|---------|--------|---------|--------|
| Build IA* Time | Baseline | ~0.5x | ~0.3-0.5x | ~0.2-0.3x |
| Forman Gradient Time | Baseline | ~Baseline | ~0.9-1.0x | ~0.9-1.0x |
| Total Time | Baseline | ~0.6x | ~0.4-0.6x | ~0.3-0.4x |
| Memory Usage | Baseline | ~Baseline | ~1.1-1.3x | ~1.1-1.3x |

---

## 9. Notes and Considerations

### Memory vs Speed Trade-off
- New code uses more memory (complete mapping) but is faster
- This is acceptable for most use cases
- Can be optimized further if needed (e.g., compressed storage)

### Parallel Build Considerations
- Parallel build should be faster for large datasets
- Sequential build may be faster for small datasets (overhead)
- Test both to find optimal choice

### Timing Accuracy
- Use high-resolution timers (`std::chrono::high_resolution_clock`)
- Run multiple times and average for consistency
- Consider system load when comparing

### File Format Considerations
- `.ply` files: Complete workflow (segmentation)
- `.off` files: Forman gradient only (performance testing)
- Ensure consistent input parameters across all runs

---

## 10. Next Steps

1. **Review this plan** and adjust as needed
2. **Start with Phase 1**: Create comparison tool and prepare git branches
3. **Implement Phase 2**: Update timing code
4. **Compile and test**: Follow the checklist
5. **Analyze results**: Create performance report

---

## 11. Appendix: File Structure

```
libTTS_public/
├── python/
│   └── tests/
│       └── cmp_ply.py          # NEW: Comparison tool
├── cpp/
│   ├── source/
│   │   ├── projects/
│   │   │   └── test_forman_gradient.cpp  # NEW: Test program
│   │   ├── iastar/
│   │   │   └── io_functions.cpp          # UPDATE: Enhanced timing
│   │   └── forman/
│   │       └── formangradient.cpp        # UPDATE: Enhanced timing
│   └── cmake-build-debug/
│       ├── xx_tts_seq_old                # OLD: Sequential
│       ├── xx_tts_pa_old                 # OLD: Parallel
│       ├── xx_tts_seq_new                # NEW: Sequential
│       └── xx_tts_pa_new                 # NEW: Parallel
└── PERFORMANCE_TESTING_PLAN.md            # This document
```

---

**Last Updated**: 2025-01-XX
**Status**: Planning Phase

