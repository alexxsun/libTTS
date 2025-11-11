#!/bin/bash
#
# Analyze test results - finds the most recent test results directory and analyzes it
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BASE_DIR="${1:-cpp/cmake-build-debug}"

echo "=========================================="
echo "Test Results Analyzer"
echo "=========================================="
echo ""

# Find most recent test results directory
LATEST_RESULTS=$(find "$BASE_DIR" -type d -name "test_results_*" 2>/dev/null | sort | tail -1)

if [ -z "$LATEST_RESULTS" ]; then
    echo "No test results directories found in $BASE_DIR"
    echo ""
    echo "To run tests:"
    echo "  ./run_tests.sh $BASE_DIR"
    exit 1
fi

echo "Found test results: $(basename $LATEST_RESULTS)"
echo ""

# Check if log files exist and have content
echo "Checking log files..."
EMPTY_COUNT=0
TOTAL_COUNT=0
HAS_CONTENT=0

for log_file in $(find "$LATEST_RESULTS" -name "*.log" -type f); do
    TOTAL_COUNT=$((TOTAL_COUNT + 1))
    if [ ! -s "$log_file" ]; then
        EMPTY_COUNT=$((EMPTY_COUNT + 1))
    else
        HAS_CONTENT=$((HAS_CONTENT + 1))
    fi
done

echo "  Total log files: $TOTAL_COUNT"
echo "  Empty files: $EMPTY_COUNT"
echo "  Files with content: $HAS_CONTENT"
echo ""

if [ $HAS_CONTENT -eq 0 ]; then
    echo "⚠ All log files are empty!"
    echo ""
    echo "Possible reasons:"
    echo "1. Executables need specific command-line arguments"
    echo "2. Tests haven't completed successfully"
    echo "3. Output redirection issue"
    echo ""
    echo "For .off files, you should use test_forman_gradient:"
    echo "  cd $BASE_DIR"
    echo "  ./test_forman_gradient <file.off> 3"
    echo ""
    echo "For .ply files, xx_tts needs -tts argument:"
    echo "  cd $BASE_DIR"
    echo "  ./xx_tts_seq_old <file.ply> <trunk_file.pts> -tts"
    echo ""
    exit 1
fi

# Run Python analysis script
echo "Running analysis..."
echo ""
python3 interpret_results.py "$LATEST_RESULTS"

