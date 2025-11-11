#!/bin/bash
# Script to compare timing between old and new versions

set -e  # Exit on error

PROJECT_DIR="/home/alex/Projects/libTTS_public"
BUILD_DIR="$PROJECT_DIR/cpp/cmake-build-debug"
TEST_FILE="$BUILD_DIR/close_stems_3_a0.01.ply"
# t109_roi_a0.010.off
TEST_FILE="$BUILD_DIR/t109_roi_a0.010.off"
TEST_LOCS="$BUILD_DIR/close_stems_3_locs.pts"

cd "$PROJECT_DIR"

echo "=== Timing Comparison Script ==="
echo ""

# Step 1: Save current changes
echo "Step 1: Stashing current changes..."
git stash push -m "New implementation without adjRelations" -- cpp/source/iastar/

# Step 2: Build old version
echo ""
echo "Step 2: Building OLD version (with adjRelations)..."
cd "$BUILD_DIR"
cmake ../ > /dev/null 2>&1
make xx_tts -j4 > /dev/null 2>&1
# rename xx_tts to xx_tts_old
mv "$BUILD_DIR/xx_tts" "$BUILD_DIR/xx_tts_old"
echo "Old version saved to: xx_tts_old"

# Step 3: Test old version
echo ""
echo "Step 3: Testing OLD version..."
if [ -f "$TEST_FILE" ] && [ -f "$TEST_LOCS" ]; then
    echo "Running old version..."
    time ./xx_tts "$TEST_FILE" "$TEST_LOCS" -tts 2>&1 | tee "$BUILD_DIR/timing_old.log"
    OLD_OUTPUT="$BUILD_DIR/close_stems_3_a0.01_lbl.ply"
    if [ -f "$OLD_OUTPUT" ]; then
        mv "$OLD_OUTPUT" "$BUILD_DIR/close_stems_3_a0.01_lbl_OLD.ply"
        echo "Old output saved to: close_stems_3_a0.01_lbl_OLD.ply"
    fi
else
    echo "WARNING: Test files not found in $BUILD_DIR"
    echo "Please ensure close_stems_3_a0.01.ply and close_stems_3_locs.pts are in the build directory"
fi

# Step 4: Restore new changes
echo ""
echo "Step 4: Restoring NEW changes..."
cd "$PROJECT_DIR"
git stash pop

# Step 5: Build new version
echo ""
echo "Step 5: Building NEW version (without adjRelations)..."
cd "$BUILD_DIR"
cmake ../ > /dev/null 2>&1
make xx_tts -j4 > /dev/null 2>&1

# Step 6: Test new version
echo ""
echo "Step 6: Testing NEW version..."
if [ -f "$TEST_FILE" ] && [ -f "$TEST_LOCS" ]; then
    echo "Running new version..."
    time ./xx_tts "$TEST_FILE" "$TEST_LOCS" -tts 2>&1 | tee "$BUILD_DIR/timing_new.log"
    NEW_OUTPUT="$BUILD_DIR/close_stems_3_a0.01_lbl.ply"
    if [ -f "$NEW_OUTPUT" ]; then
        mv "$NEW_OUTPUT" "$BUILD_DIR/close_stems_3_a0.01_lbl_NEW.ply"
        echo "New output saved to: close_stems_3_a0.01_lbl_NEW.ply"
    fi
else
    echo "WARNING: Test files not found in $BUILD_DIR"
fi

# Step 7: Compare results
echo ""
echo "=== Timing Comparison ==="
echo ""
echo "OLD VERSION (with adjRelations):"
if [ -f "$BUILD_DIR/timing_old.log" ]; then
    echo "  Data structure building:"
    grep -E "(build IA\*|Complete mapping storage|adj phase)" "$BUILD_DIR/timing_old.log" | head -5
    echo "  storeFullStar:"
    grep -E "Tops computed" "$BUILD_DIR/timing_old.log" || echo "    Not found"
    echo "  Overall times:"
    grep -E "(grow time|label time)" "$BUILD_DIR/timing_old.log" || echo "    Not found"
fi

echo ""
echo "NEW VERSION (without adjRelations):"
if [ -f "$BUILD_DIR/timing_new.log" ]; then
    echo "  Data structure building:"
    grep -E "(build IA\*|Complete mapping storage|adj phase)" "$BUILD_DIR/timing_new.log" | head -5
    echo "  storeFullStar:"
    grep -E "Tops computed" "$BUILD_DIR/timing_new.log" || echo "    Not found"
    echo "  Overall times:"
    grep -E "(grow time|label time)" "$BUILD_DIR/timing_new.log" || echo "    Not found"
fi

echo ""
echo "=== Summary ==="
echo "Old timing log: $BUILD_DIR/timing_old.log"
echo "New timing log: $BUILD_DIR/timing_new.log"
echo ""
echo "To see full comparison, run:"
echo "  diff $BUILD_DIR/timing_old.log $BUILD_DIR/timing_new.log"
echo ""
echo "Done!"

