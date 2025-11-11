#!/bin/bash
#
# Script to compile 4 variants of xx_tts:
# - xx_tts_seq_old: Old code (main branch), sequential build
# - xx_tts_pa_old: Old code (main branch), sequential build (same as seq_old)
# - xx_tts_seq_new: New code (improve_iastar branch), sequential build
# - xx_tts_pa_new: New code (improve_iastar branch), parallel build
#
# Usage:
#   ./compile_variants.sh [build_dir]
#
# Default build_dir: cpp/cmake-build-debug
#
# Prerequisites:
#   - main branch: old code with adjRelations, enhanced timing
#   - improve_iastar branch: new code with completeCoboundaryTop, enhanced timing
#

set -e  # Exit on error

BUILD_DIR="${1:-cpp/cmake-build-debug}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "Compiling 4 xx_tts Variants"
echo "=========================================="
echo "Build directory: $BUILD_DIR"
echo ""

# Check if git is available
if ! command -v git &> /dev/null; then
    echo "Error: git is not installed"
    exit 1
fi

# Get current branch
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
echo "Current branch: $CURRENT_BRANCH"
echo ""

# Check if required branches exist
if ! git show-ref --verify --quiet refs/heads/main; then
    echo "Error: main branch not found"
    exit 1
fi

if ! git show-ref --verify --quiet refs/heads/improve_iastar; then
    echo "Error: improve_iastar branch not found"
    exit 1
fi

# Function to compile a variant
compile_variant() {
    local variant_name=$1
    local use_parallel=$2
    
    echo "----------------------------------------"
    echo "Compiling: $variant_name"
    echo "  Parallel: $use_parallel"
    echo "----------------------------------------"
    
    # Create build directory
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    # Clean previous build
    rm -rf CMakeCache.txt CMakeFiles/
    
    # Configure CMake
    if [ "$use_parallel" = "true" ]; then
        cmake ../ -DCMAKE_BUILD_TYPE=Release -DUSE_PARALLEL_BUILD=ON
    else
        cmake ../ -DCMAKE_BUILD_TYPE=Release -DUSE_PARALLEL_BUILD=OFF
    fi
    
    # Build
    make xx_tts -j$(nproc) || {
        echo "✗ Build failed for $variant_name"
        cd "$SCRIPT_DIR"
        return 1
    }
    
    # Rename executable
    if [ -f xx_tts ]; then
        mv xx_tts "$variant_name"
        echo "✓ Created $variant_name"
    else
        echo "✗ Executable not found after build"
        cd "$SCRIPT_DIR"
        return 1
    fi
    
    cd "$SCRIPT_DIR"
    echo ""
}

# Compile old code variants from main branch
echo "=========================================="
echo "Compiling OLD code variants (main branch)"
echo "=========================================="
echo ""

# Switch to main branch
git checkout main
echo "Switched to main branch"
echo ""

# Old code only supports sequential (buildDataStructure() only)
compile_variant "xx_tts_seq_old" "false"
# pa_old is the same as seq_old for old code
cp "$BUILD_DIR/xx_tts_seq_old" "$BUILD_DIR/xx_tts_pa_old" 2>/dev/null || true
if [ -f "$BUILD_DIR/xx_tts_pa_old" ]; then
    echo "✓ Created xx_tts_pa_old (same as seq_old for old code)"
    echo ""
fi

# Compile new code variants from improve_iastar branch
echo "=========================================="
echo "Compiling NEW code variants (improve_iastar branch)"
echo "=========================================="
echo ""

# Switch to improve_iastar branch
git checkout improve_iastar
echo "Switched to improve_iastar branch"
echo ""

compile_variant "xx_tts_seq_new" "false"
compile_variant "xx_tts_pa_new" "true"

# Switch back to original branch
git checkout "$CURRENT_BRANCH"
echo ""
echo "Switched back to $CURRENT_BRANCH branch"
echo ""

echo "=========================================="
echo "Compilation Summary"
echo "=========================================="
cd "$BUILD_DIR"
echo ""
echo "Created executables:"
ls -lh xx_tts_* 2>/dev/null || echo "No executables found"
echo ""
echo "Note: xx_tts_pa_old is the same as xx_tts_seq_old"
echo "      (old code only supports sequential build)"
echo ""
echo "Done!"

