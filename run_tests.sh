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
TRUNK_FILE="$TEST_DIR/close_stems_3_locs.pts"
OFF_FILES=(
    "$TEST_DIR/close_stems_3_a0.010.off"
    "$TEST_DIR/aoi_thin_low_pts_a0.010.off"
    "$TEST_DIR/tree_228_veg_xyz_as_0.01.off"
    "$TEST_DIR/t109_roi_a0.010.off"
)

# Executable variants
XX_TTS_VARIANTS=(
    "xx_tts_seq_old"
    "xx_tts_pa_old"
    "xx_tts_seq_new"
    "xx_tts_pa_new"
)

TEST_FORMAN_VARIANTS=(
    "test_forman_gradient_seq_old"
    "test_forman_gradient_pa_old"
    "test_forman_gradient_seq_new"
    "test_forman_gradient_pa_new"
)

# Check executables
echo "Checking executables..."
echo "  xx_tts variants:"
for exe in "${XX_TTS_VARIANTS[@]}"; do
    if [ -f "$TEST_DIR/$exe" ]; then
        echo "    ✓ $exe"
    else
        echo "    ✗ $exe (missing)"
    fi
done

echo "  test_forman_gradient variants:"
for exe in "${TEST_FORMAN_VARIANTS[@]}"; do
    if [ -f "$TEST_DIR/$exe" ]; then
        echo "    ✓ $exe"
    else
        echo "    ✗ $exe (missing)"
    fi
done
echo ""

# Create output directory
OUTPUT_DIR="$TEST_DIR/test_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Test .ply file (complete workflow with xx_tts)
if [ -f "$PLY_FILE" ] && [ -f "$TRUNK_FILE" ]; then
    echo "=========================================="
    echo "Testing .ply file (Complete Workflow)"
    echo "=========================================="
    echo ""
    
    PLY_OUTPUT_DIR="$OUTPUT_DIR/ply_tests"
    mkdir -p "$PLY_OUTPUT_DIR"
    
    cd "$TEST_DIR"
    for exe in "${XX_TTS_VARIANTS[@]}"; do
        if [ ! -f "$exe" ]; then
            echo "⚠ Skipping $exe (not found)"
            continue
        fi
        
        # Use absolute path for logfile since we're in TEST_DIR
        logfile="$(cd "$SCRIPT_DIR" && pwd)/$PLY_OUTPUT_DIR/${exe}_output.log"
        mkdir -p "$(dirname "$logfile")"
        
        echo "Running: $exe on $(basename $PLY_FILE)..."
        "./$exe" "$(basename "$PLY_FILE")" "$(basename "$TRUNK_FILE")" -tts 2>&1 | tee "$logfile"
        echo ""
    done
    cd "$SCRIPT_DIR"
    
    echo "To compare segmentation results:"
    echo "  python python/tests/cmp_ply.py <output1.ply> <output2.ply>"
    echo ""
else
    echo "=========================================="
    echo "Skipping .ply file test"
    echo "  PLY file: $([ -f "$PLY_FILE" ] && echo "✓" || echo "✗") $(basename "$PLY_FILE")"
    echo "  Trunk file: $([ -f "$TRUNK_FILE" ] && echo "✓" || echo "✗") $(basename "$TRUNK_FILE")"
    echo "=========================================="
    echo ""
fi

# Test .off files (Forman gradient with test_forman_gradient)
echo "=========================================="
echo "Testing .off files (Forman Gradient)"
echo "=========================================="
echo ""

OFF_OUTPUT_DIR="$OUTPUT_DIR/off_tests"
mkdir -p "$OFF_OUTPUT_DIR"

cd "$TEST_DIR"
for off_file in "${OFF_FILES[@]}"; do
    off_basename="$(basename "$off_file")"
    if [ ! -f "$off_basename" ]; then
        echo "⚠ Skipping $off_basename (not found)"
        continue
    fi
    
    echo "----------------------------------------"
    echo "File: $off_basename"
    echo "----------------------------------------"
    
    FILE_OUTPUT_DIR="$OFF_OUTPUT_DIR/$(basename "$off_basename" .off)"
    mkdir -p "$FILE_OUTPUT_DIR"
    
    for exe in "${TEST_FORMAN_VARIANTS[@]}"; do
        if [ ! -f "$exe" ]; then
            echo "⚠ Skipping $exe (not found)"
            continue
        fi
        
        # Use absolute path for logfile since we're in TEST_DIR (same pattern as ply_tests)
        logfile="$(cd "$SCRIPT_DIR" && pwd)/$FILE_OUTPUT_DIR/${exe}_output.log"
        mkdir -p "$(dirname "$logfile")"
        
        echo "Running: $exe on $off_basename..."
        "./$exe" "$off_basename" 3 2>&1 | tee "$logfile"
        echo ""
    done
    echo ""
done
cd "$SCRIPT_DIR"

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
