# Testing Workflow Documentation

## Overview

This document describes the complete workflow for compiling variants, running performance tests, and generating performance reports. The workflow consists of three main scripts that automate the entire testing process.

---

## Workflow Steps

### Step 1: Compile Variants (`compile_variants.sh`)

**Purpose**: Compiles all executable variants needed for performance comparison.

**What it does**:
- Compiles 4 variants of `xx_tts` (old/new code, sequential/parallel)
- Compiles 4 variants of `test_forman_gradient` (old/new code, sequential/parallel)
- Automatically switches between `main` and `improve_iastar` git branches
- Handles CMake configuration for sequential vs parallel builds

**Usage**:
```bash
./compile_variants.sh [build_dir]
```

**Default**: `cpp/cmake-build-debug`

**Prerequisites**:
- Git repository with `main` and `improve_iastar` branches
- CMake and build tools installed
- OpenMP support (for parallel builds)

**Output**:
- `xx_tts_seq_old` - Old code (main branch), sequential build
- `xx_tts_pa_old` - Old code (main branch), sequential build (same as seq_old)
- `xx_tts_seq_new` - New code (improve_iastar branch), sequential build
- `xx_tts_pa_new` - New code (improve_iastar branch), parallel build
- `test_forman_gradient_seq_old` - Old code (main branch), sequential build
- `test_forman_gradient_pa_old` - Old code (main branch), sequential build
- `test_forman_gradient_seq_new` - New code (improve_iastar branch), sequential build
- `test_forman_gradient_pa_new` - New code (improve_iastar branch), parallel build

**How it works**:
1. Checks for required git branches (`main` and `improve_iastar`)
2. For each variant:
   - Switches to appropriate branch
   - Configures CMake with appropriate flags (`USE_PARALLEL_BUILD=ON/OFF`)
   - Compiles the target executable
   - Renames executable with variant suffix
3. Returns to original branch

**Example**:
```bash
cd /home/alex/Projects/libTTS_public
./compile_variants.sh cpp/cmake-build-debug
```

---

### Step 2: Run Tests (`run_tests.sh`)

**Purpose**: Runs performance tests on all test files and collects timing data.

**What it does**:
- Tests `.ply` files with complete workflow (using `xx_tts` variants)
- Tests `.off` files with Forman gradient only (using `test_forman_gradient` variants)
- Collects timing data in log files
- Organizes output by test file and variant

**Usage**:
```bash
./run_tests.sh [test_dir]
```

**Default**: `cpp/cmake-build-debug`

**Prerequisites**:
- All executable variants compiled (from Step 1)
- Test files in the test directory:
  - `close_stems_3_a0.01.ply` (with `close_stems_3_locs.pts` trunk file)
  - `close_stems_3_a0.010.off`
  - `aoi_thin_low_pts_a0.010.off`
  - `tree_228_veg_xyz_as_0.01.off`
  - `t109_roi_a0.010.off`

**Output Structure**:
```
test_results_YYYYMMDD_HHMMSS/
├── ply_tests/
│   ├── xx_tts_seq_old_output.log
│   ├── xx_tts_pa_old_output.log
│   ├── xx_tts_seq_new_output.log
│   └── xx_tts_pa_new_output.log
└── off_tests/
    ├── close_stems_3_a0.010/
    │   ├── test_forman_gradient_seq_old.log
    │   ├── test_forman_gradient_pa_old.log
    │   ├── test_forman_gradient_seq_new.log
    │   └── test_forman_gradient_pa_new.log
    ├── aoi_thin_low_pts_a0.010/
    │   └── ...
    └── ...
```

**How it works**:
1. Checks for required executables
2. Creates timestamped output directory
3. For `.ply` files:
   - Runs each `xx_tts` variant with `-tts` flag and trunk file
   - Captures output to log files using `tee`
4. For `.off` files:
   - Runs each `test_forman_gradient` variant
   - Organizes logs by test file name
   - Captures output to log files using `tee`
5. All output is saved with both stdout and stderr captured

**Example**:
```bash
cd /home/alex/Projects/libTTS_public
./run_tests.sh cpp/cmake-build-debug
```

**Log File Format**:
Each log file contains:
- File statistics (vertices, top simplexes)
- Detailed timing breakdown:
  - File I/O time
  - Vertex reading time
  - Cell reading time
  - Build IA* time (sequential/parallel)
  - Total reading time
  - Forman gradient computation time
  - Total execution time
- Memory usage (KB and MB)

---

### Step 3: Generate Performance Report (`generate_performance_report.py`)

**Purpose**: Parses test logs and generates comprehensive performance reports with figures.

**What it does**:
- Parses timing data from all log files
- Identifies old vs new and sequential vs parallel variants
- Generates comparison plots using matplotlib
- Creates text summary report

**Usage**:
```bash
python generate_performance_report.py <test_results_dir>
```

**Example**:
```bash
python generate_performance_report.py cpp/cmake-build-debug/test_results_20251111_123424
```

**Prerequisites**:
- Python 3 with matplotlib and numpy (optional, for plots)
- Test results directory from Step 2

**Output Files** (in `test_results_dir/report/`):
- `sequential_comparison.png` - Old vs new sequential comparison (6 subplots)
- `parallel_comparison.png` - Old vs new parallel comparison (6 subplots)
- `parallel_vs_sequential.png` - New code sequential vs parallel (6 subplots)
- `memory_usage.png` - Memory usage comparison (4 subplots)
- `performance_report.txt` - Text summary with all metrics

**How it works**:
1. Scans test results directory for log files
2. Parses timing data using regex patterns:
   - File I/O, vertex reading, cell reading
   - Build IA* time (sequential/parallel)
   - Total reading time
   - Forman gradient computation time
   - Total execution time
   - Memory usage
3. Identifies variants:
   - Old: `*_old` suffix
   - New: `*_new` suffix
   - Sequential: `*_seq_*` or sequential build flag
   - Parallel: `*_pa_*` or parallel build flag
4. Generates comparison plots:
   - Time comparisons (bar charts)
   - Speedup analysis (old/new ratios)
   - Scalability analysis (time vs file size)
   - Memory usage comparison
5. Creates text report with summary tables

**Plot Details**:

**`sequential_comparison.png`** (Old vs New Sequential):
- Subplot 1: Build IA* time
- Subplot 2: Total reading time
- Subplot 3: Top complex generation time
- Subplot 4: Forman gradient calculation time
- Subplot 5: Detailed speedup breakdown (Build IA, Total Reading, Forman Gradient, Total Time)
- Subplot 6: Scalability (time vs # vertices, # faces, avg faces per vertex)

**`parallel_comparison.png`** (Old vs New Parallel):
- Same structure as sequential comparison

**`parallel_vs_sequential.png`** (New Code Sequential vs Parallel):
- Subplot 1: Build IA* time
- Subplot 2: Total reading time
- Subplot 3: Forman gradient calculation time
- Subplot 4: Total time
- Subplot 5: Parallel speedup breakdown (Build IA, Total Reading, Forman Gradient, Total Time)
- Subplot 6: Scalability comparison

**`memory_usage.png`** (Memory Usage Comparison):
- Subplot 1: Sequential old vs new
- Subplot 2: Parallel old vs new
- Subplot 3: New code sequential vs parallel
- Subplot 4: Overall comparison

---

## Complete Workflow Example

```bash
# Step 1: Compile all variants
cd /home/alex/Projects/libTTS_public
./compile_variants.sh cpp/cmake-build-debug

# Step 2: Run tests
./run_tests.sh cpp/cmake-build-debug

# Step 3: Generate report (use the timestamped directory from Step 2)
python generate_performance_report.py cpp/cmake-build-debug/test_results_20251111_123424

# View results
ls -lh cpp/cmake-build-debug/test_results_*/report/
```

---

## Script Details

### `compile_variants.sh`

**Key Features**:
- Automatic git branch switching
- CMake configuration for parallel/sequential builds
- Executable renaming with variant suffixes
- Error handling and validation

**Configuration**:
- Build type: `Release` (hardcoded)
- Parallel flag: `USE_PARALLEL_BUILD=ON/OFF` (controlled by variant)
- Build directory: Configurable via argument

**Dependencies**:
- Git
- CMake
- C++ compiler with OpenMP support

---

### `run_tests.sh`

**Key Features**:
- Automatic test file detection
- Organized output directory structure
- Captures both stdout and stderr using `tee`
- Handles missing executables gracefully

**Test Files**:
- `.ply` files: Complete workflow (segmentation)
- `.off` files: Forman gradient only

**Command-line Arguments**:
- `xx_tts` variants: `./xx_tts <ply_file> <trunk_file> -tts`
- `test_forman_gradient` variants: `./test_forman_gradient <off_file>`

**Output Organization**:
- Separate directories for `.ply` and `.off` tests
- Log files named by variant and test file
- Timestamped output directory

---

### `generate_performance_report.py`

**Key Features**:
- Automatic log file parsing
- Variant identification (old/new, seq/parallel)
- Comprehensive plot generation
- Text summary report

**Parsing Patterns**:
- Timing values: `X.XXX s` format
- Memory usage: `Memory usage (KB): XXXX` and `Memory usage (MB): X.XX`
- File statistics: `Vertices: X`, `Total top simplexes: X`

**Plot Generation**:
- Uses matplotlib with non-interactive backend (`Agg`)
- Creates publication-quality figures
- Handles missing data gracefully
- Log-log axes for scalability plots

**Error Handling**:
- Checks for matplotlib availability
- Falls back to text-only report if matplotlib not available
- Validates input directory existence
- Handles missing or empty log files

---

## Troubleshooting

### Compilation Issues

**Problem**: `compile_variants.sh` fails to find branches
- **Solution**: Ensure `main` and `improve_iastar` branches exist
- **Check**: `git branch -a`

**Problem**: CMake configuration fails
- **Solution**: Check CMake and compiler installation
- **Check**: `cmake --version`, `g++ --version`

**Problem**: OpenMP not found
- **Solution**: Install OpenMP development packages
- **Ubuntu/Debian**: `sudo apt-get install libomp-dev`
- **Arch/Manjaro**: `sudo pacman -S openmp`

### Test Execution Issues

**Problem**: Executables not found
- **Solution**: Run `compile_variants.sh` first
- **Check**: `ls -lh cpp/cmake-build-debug/*_seq_* cpp/cmake-build-debug/*_pa_*`

**Problem**: Test files not found
- **Solution**: Ensure test files are in the test directory
- **Check**: `ls -lh cpp/cmake-build-debug/*.ply cpp/cmake-build-debug/*.off`

**Problem**: Empty log files
- **Solution**: Check that executables produce output
- **Check**: Run executables manually to verify output

### Report Generation Issues

**Problem**: matplotlib not found
- **Solution**: Install matplotlib and numpy
- **Command**: `python3 -m pip install --user matplotlib numpy`
- **Note**: Script will generate text-only report if matplotlib unavailable

**Problem**: No results found
- **Solution**: Check that log files contain timing data
- **Check**: `grep "build IA" cpp/cmake-build-debug/test_results_*/*/*.log`

**Problem**: Missing data in plots
- **Solution**: Check that log files are properly formatted
- **Check**: Verify timing output format matches expected patterns

---

## Tips and Best Practices

1. **Always compile variants first**: Run `compile_variants.sh` before `run_tests.sh`
2. **Check prerequisites**: Ensure all test files and executables exist
3. **Review log files**: Check log files for errors before generating reports
4. **Use absolute paths**: Scripts handle paths automatically, but be aware of working directory
5. **Keep git branches clean**: Ensure `main` and `improve_iastar` branches are up-to-date
6. **Monitor disk space**: Test results can be large, especially with many test files

---

## Output Interpretation

### Timing Metrics

- **Build IA***: Time to build the incidence array data structure
- **Total Reading**: Time for complete file reading and data structure building
- **Forman Gradient**: Time for Forman gradient computation
- **Total Time**: Complete execution time

### Speedup Calculation

Speedup = `old_time / new_time`
- Speedup > 1.0: New code is faster
- Speedup < 1.0: New code is slower (shouldn't happen)

### Scalability Analysis

- **Time vs # Vertices**: How performance scales with mesh size
- **Time vs # Faces**: How performance scales with complexity
- **Time vs Avg Faces per Vertex**: How performance scales with connectivity

---

**Last Updated**: 2025-Nov-13  
**Status**: ✅ Complete and tested

