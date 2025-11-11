#!/bin/bash
#
# Helper script to collect timing data from test outputs
# Run this after running tests to extract timing information
#

set -e

BASE_DIR="${1:-cpp/cmake-build-debug}"

echo "=========================================="
echo "Timing Data Collection Helper"
echo "=========================================="
echo ""

echo "This script helps you collect timing data from test outputs."
echo ""
echo "To collect data manually:"
echo ""
echo "1. For each test file, run:"
echo "   ./test_forman_gradient <file.off> 3 > <file>_<variant>.log"
echo ""
echo "2. Extract key metrics from logs:"
echo "   - Vertices: Look for 'Vertices: X'"
echo "   - Top Simplexes: Look for 'Total top simplexes: X'"
echo "   - Build IA*: Look for 'build IA* (sequential/parallel): X.XXX s'"
echo "   - Forman Gradient: Look for 'Total Forman gradient time: X.XXX s'"
echo "   - Total Time: Look for 'Total time: X.XXX s'"
echo ""
echo "3. Fill in the PERFORMANCE_REPORT.md template"
echo ""
echo "Example commands to run tests:"
echo "----------------------------------------"
echo ""

# Check if executables exist
if [ -f "$BASE_DIR/test_forman_gradient" ]; then
    echo "Found test_forman_gradient program"
    echo ""
    echo "Run tests like this:"
    echo ""
    
    for file in "$BASE_DIR"/*.off; do
        if [ -f "$file" ]; then
            filename=$(basename "$file")
            echo "  # $filename"
            echo "  cd $BASE_DIR"
            echo "  ./test_forman_gradient $filename 3 > ${filename%.off}_test.log"
            echo ""
        fi
    done
else
    echo "test_forman_gradient not found in $BASE_DIR"
    echo "Make sure you've compiled it first"
fi

echo ""
echo "Or use the xx_tts executables (but they run full workflow):"
echo ""

for variant in seq_old pa_old seq_new pa_new; do
    if [ -f "$BASE_DIR/xx_tts_$variant" ]; then
        echo "  Found: xx_tts_$variant"
    fi
done

echo ""
echo "=========================================="
echo "Data Collection Complete!"
echo "=========================================="

