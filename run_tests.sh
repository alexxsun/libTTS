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

# Also check for test_forman_gradient variants
echo ""
echo "Checking test_forman_gradient executables..."
TEST_EXECUTABLES=(
    "test_forman_gradient_seq_old"
    "test_forman_gradient_pa_old"
    "test_forman_gradient_seq_new"
    "test_forman_gradient_pa_new"
)

for test_exe in "${TEST_EXECUTABLES[@]}"; do
    if [ -f "$TEST_DIR/$test_exe" ]; then
        echo "  ✓ Found $test_exe"
    else
        echo "  ✗ Missing $test_exe"
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
    
    # Ensure output directory exists and convert to absolute path
    mkdir -p "$(dirname "$output_file")"
    local abs_output_file
    if [[ "$output_file" = /* ]]; then
        abs_output_file="$output_file"
    else
        # Convert relative path to absolute
        abs_output_file="$(cd "$(dirname "$output_file")" && pwd)/$(basename "$output_file")"
    fi
    
    # Run test (redirect both stdout and stderr using tee)
    cd "$TEST_DIR"
    "./$exe" "$input_file" 2>&1 | tee "$abs_output_file" || true
    cd "$SCRIPT_DIR"
    
    if [ -s "$abs_output_file" ]; then
        echo "  ✓ Output saved to $abs_output_file ($(wc -l < "$abs_output_file") lines)"
    else
        echo "  ⚠ Output file is empty - executable may need different arguments"
    fi
}

# Test .ply file (complete workflow)
if [ -f "$PLY_FILE" ]; then
    echo "=========================================="
    echo "Testing .ply file (Complete Workflow)"
    echo "=========================================="
    echo ""
    
    PLY_OUTPUT_DIR="$OUTPUT_DIR/ply_tests"
    mkdir -p "$PLY_OUTPUT_DIR"
    
    # Check for trunk file (required for -tts)
    TRUNK_FILE="$TEST_DIR/close_stems_3_locs.pts"
    if [ ! -f "$TRUNK_FILE" ]; then
        echo "⚠ Warning: Trunk file not found: $TRUNK_FILE"
        echo "  xx_tts requires -tts argument with trunk file for .ply files"
        echo "  Skipping .ply file tests or they may produce empty output"
        echo ""
    fi
    
    for exe in "${EXECUTABLES[@]}"; do
        output_file="$PLY_OUTPUT_DIR/${exe}_output.log"
        # Ensure output directory exists and convert to absolute path
        mkdir -p "$(dirname "$output_file")"
        abs_output_file=""
        if [[ "$output_file" = /* ]]; then
            abs_output_file="$output_file"
        else
            abs_output_file="$(cd "$(dirname "$output_file")" && pwd)/$(basename "$output_file")"
        fi
        
        echo "Running: $exe on $(basename $PLY_FILE)..."
        if [ -f "$TEST_DIR/$exe" ] && [ -f "$TRUNK_FILE" ]; then
            cd "$TEST_DIR"
            "./$exe" "$PLY_FILE" "$TRUNK_FILE" -tts 2>&1 | tee "$abs_output_file" || true
            cd "$SCRIPT_DIR"
            if [ -s "$abs_output_file" ]; then
                echo "  ✓ Output saved to $abs_output_file ($(wc -l < "$abs_output_file") lines)"
            else
                echo "  ⚠ Output file is empty"
            fi
        else
            echo "  ⚠ Skipping (executable or trunk file not found)"
        fi
        echo ""
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
    
    # Test with test_forman_gradient variants
    # Map xx_tts variants to test_forman_gradient variants
    declare -A TEST_EXECUTABLES=(
        ["xx_tts_seq_old"]="test_forman_gradient_seq_old"
        ["xx_tts_pa_old"]="test_forman_gradient_pa_old"
        ["xx_tts_seq_new"]="test_forman_gradient_seq_new"
        ["xx_tts_pa_new"]="test_forman_gradient_pa_new"
    )
    
    echo "Using test_forman_gradient variants..."
    for exe in "${EXECUTABLES[@]}"; do
        test_exe="${TEST_EXECUTABLES[$exe]}"
        if [ -z "$test_exe" ]; then
            echo "  ⚠ No test_forman_gradient variant for $exe, skipping"
            continue
        fi
        
        output_file="$FILE_OUTPUT_DIR/${test_exe}_output.log"
        # Ensure output directory exists and convert to absolute path
        mkdir -p "$(dirname "$output_file")"
        if [[ "$output_file" = /* ]]; then
            abs_output_file="$output_file"
        else
            abs_output_file="$(cd "$(dirname "$output_file")" && pwd)/$(basename "$output_file")"
        fi
        
        echo "Running: $test_exe on $(basename $off_file)..."
        if [ -f "$TEST_DIR/$test_exe" ]; then
            cd "$TEST_DIR"
            "./$test_exe" "$off_file" 3 2>&1 | tee "$abs_output_file" || true
            cd "$SCRIPT_DIR"
            if [ -s "$abs_output_file" ]; then
                echo "  ✓ Output saved to $abs_output_file ($(wc -l < "$abs_output_file") lines)"
            else
                echo "  ⚠ Output file is empty"
            fi
        else
            echo "  ⚠ Executable not found: $test_exe"
            echo "     Run ./compile_variants.sh to build all variants"
        fi
        echo ""
    done
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

