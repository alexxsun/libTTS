#!/bin/bash
#
# Complete performance analysis workflow
# This script helps complete all remaining analysis steps
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BASE_DIR="${1:-cpp/cmake-build-debug}"

echo "=========================================="
echo "Performance Analysis - Complete Workflow"
echo "=========================================="
echo ""

# Step 1: Compare segmentation results
echo "Step 1: Comparing segmentation results..."
echo "----------------------------------------"

PLY_OLD="$BASE_DIR/close_stems_3_a0.01_lbl_OLD.ply"
PLY_NEW="$BASE_DIR/close_stems_3_a0.01_lbl_NEW.ply"

if [ -f "$PLY_OLD" ] && [ -f "$PLY_NEW" ]; then
    echo "Found segmentation files:"
    echo "  Old: $PLY_OLD"
    echo "  New: $PLY_NEW"
    echo ""
    
    # Try to activate venv and run comparison
    if [ -f ~/xx_pyvenvs/treemapping_project/bin/activate.fish ]; then
        echo "Running comparison (using fish venv)..."
        fish -c "source ~/xx_pyvenvs/treemapping_project/bin/activate.fish && python3 python/tests/cmp_ply.py '$PLY_OLD' '$PLY_NEW'" || {
            echo "Note: Could not run comparison (venv issue)."
            echo "      You can run manually:"
            echo "      . ~/xx_pyvenvs/treemapping_project/bin/activate.fish"
            echo "      python python/tests/cmp_ply.py $PLY_OLD $PLY_NEW"
        }
    else
        echo "Note: Virtual environment not found. Run comparison manually:"
        echo "      . ~/xx_pyvenvs/treemapping_project/bin/activate.fish"
        echo "      python python/tests/cmp_ply.py $PLY_OLD $PLY_NEW"
    fi
else
    echo "Segmentation files not found. Skipping comparison."
fi

echo ""
echo "Step 2: Collecting timing data..."
echo "----------------------------------------"

# Run performance report
python3 create_performance_report.py "$BASE_DIR"

echo ""
echo "Step 3: Next steps..."
echo "----------------------------------------"
echo ""
echo "To complete the analysis:"
echo ""
echo "1. Run tests with enhanced timing (if not done):"
echo "   ./run_tests.sh $BASE_DIR"
echo ""
echo "2. Compare segmentation results:"
echo "   . ~/xx_pyvenvs/treemapping_project/bin/activate.fish"
echo "   python python/tests/cmp_ply.py $BASE_DIR/close_stems_3_a0.01_lbl_OLD.ply $BASE_DIR/close_stems_3_a0.01_lbl_NEW.ply"
echo ""
echo "3. For each test file, collect:"
echo "   - Vertex count"
echo "   - Top simplex count"
echo "   - Build IA* time (old vs new, seq vs parallel)"
echo "   - Forman gradient time"
echo "   - Total time"
echo ""
echo "4. Create scalability analysis:"
echo "   - Time vs # vertices"
echo "   - Time per vertex metrics"
echo "   - Speedup calculations (old vs new)"
echo ""
echo "5. Generate final report with:"
echo "   - Performance comparison tables"
echo "   - Scalability analysis"
echo "   - Correctness verification (IoU results)"
echo ""

