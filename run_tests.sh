#!/bin/bash
#
# Script to run performance tests on all test files
# Collects timing data and generates summary reports
#
# Usage:
#   ./run_tests.sh [test_dir]
#
# Default test_dir: cpp/cmake-build-debug
#

set -e

TEST_DIR="${1:-cpp/cmake-build-debug}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "Performance Testing Suite"
echo "=========================================="
echo "Test directory: $TEST_DIR"
echo ""

# Test files
PLY_FILE="$TEST_DIR/close_stems_3_a0.01.ply"
OFF_FILES=(
    "$TEST_DIR/close_stems_3_a0.010.off"
    "$TEST_DIR/aoi_thin_low_pts_a0.010.off"
    "$TEST_DIR/tree_228_veg_xyz_as_0.01.off"
)

# Executables to test
EXECUTABLES=(
    "xx_tts_seq_old"
    "xx_tts_pa_old"
    "xx_tts_seq_new"
    "xx_tts_pa_new"
)

# Check if executables exist
echo "Checking executables..."
for exe in "${EXECUTABLES[@]}"; do
    if [ -f "$TEST_DIR/$exe" ]; then
        echo "  ✓ Found $exe"
    else
        echo "  ✗ Missing $exe"
    fi
done
echo ""

# Create output directory
OUTPUT_DIR="$TEST_DIR/test_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Function to run test and capture output
run_test() {
    local exe=$1
    local input_file=$2
    local output_file=$3
    
    echo "Running: $exe on $(basename $input_file)..."
    
    if [ ! -f "$TEST_DIR/$exe" ]; then
        echo "  ✗ Executable not found, skipping"
        return 1
    fi
    
    if [ ! -f "$input_file" ]; then
        echo "  ✗ Input file not found, skipping"
        return 1
    fi
    
    # Run test (redirect both stdout and stderr)
    cd "$TEST_DIR"
    "./$exe" "$input_file" > "$output_file" 2>&1 || true
    cd "$SCRIPT_DIR"
    
    echo "  ✓ Output saved to $(basename $output_file)"
}

# Test .ply file (complete workflow)
if [ -f "$PLY_FILE" ]; then
    echo "=========================================="
    echo "Testing .ply file (Complete Workflow)"
    echo "=========================================="
    echo ""
    
    PLY_OUTPUT_DIR="$OUTPUT_DIR/ply_tests"
    mkdir -p "$PLY_OUTPUT_DIR"
    
    for exe in "${EXECUTABLES[@]}"; do
        output_file="$PLY_OUTPUT_DIR/${exe}_output.log"
        run_test "$exe" "$PLY_FILE" "$output_file"
        
        # Extract output .ply file if it exists
        # (This depends on the actual workflow output)
    done
    
    echo ""
    echo "To compare segmentation results:"
    echo "  . ~/xx_pyvenvs/treemapping_project/bin/activate.fish"
    echo "  python python/tests/cmp_ply.py <output1.ply> <output2.ply>"
    echo ""
else
    echo "=========================================="
    echo "Skipping .ply file test (file not found)"
    echo "=========================================="
    echo ""
fi

# Test .off files (Forman gradient only)
echo "=========================================="
echo "Testing .off files (Forman Gradient)"
echo "=========================================="
echo ""

OFF_OUTPUT_DIR="$OUTPUT_DIR/off_tests"
mkdir -p "$OFF_OUTPUT_DIR"

for off_file in "${OFF_FILES[@]}"; do
    if [ ! -f "$off_file" ]; then
        echo "Skipping $(basename $off_file) (not found)"
        continue
    fi
    
    echo "----------------------------------------"
    echo "File: $(basename $off_file)"
    echo "----------------------------------------"
    
    FILE_OUTPUT_DIR="$OFF_OUTPUT_DIR/$(basename $off_file .off)"
    mkdir -p "$FILE_OUTPUT_DIR"
    
    # Test with test_forman_gradient if available
    if [ -f "$TEST_DIR/test_forman_gradient" ]; then
        echo "Using test_forman_gradient program..."
        for exe in "${EXECUTABLES[@]}"; do
            # For test_forman_gradient, we need to compile variants
            # For now, just test with the main test program
            output_file="$FILE_OUTPUT_DIR/test_forman_gradient.log"
            cd "$TEST_DIR"
            "./test_forman_gradient" "$off_file" 3 > "$output_file" 2>&1 || true
            cd "$SCRIPT_DIR"
        done
    else
        echo "test_forman_gradient not found, using xx_tts executables..."
        for exe in "${EXECUTABLES[@]}"; do
            output_file="$FILE_OUTPUT_DIR/${exe}_output.log"
            run_test "$exe" "$off_file" "$output_file"
        done
    fi
    echo ""
done

echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Next steps:"
echo "1. Review timing data in output files"
echo "2. Compare segmentation results for .ply file"
echo "3. Analyze scalability (file size vs performance)"
echo "4. Generate performance report"
echo ""

