# Git Branch Setup Plan for Performance Testing

## Current Situation

- **`main` branch**: Has old code with `adjRelations`, but no standardized timing output
- **`improve_iastar` branch**: Has new code with `completeCoboundaryTop` (uncommitted changes), with enhanced timing
- Both branches are at the same commit (`f31c77c`)
- All differences are currently uncommitted in `improve_iastar`

## Goal

Set up two branches for performance comparison:
1. **Old code branch** (`main`): Keep `adjRelations`, add timing, use `buildDataStructure()` only
2. **New code branch** (`improve_iastar`): Keep `completeCoboundaryTop`, keep timing, support both sequential and parallel

## Step-by-Step Plan

### Phase 1: Save Current Work in `improve_iastar`

1. **Commit all changes in `improve_iastar`**:
   ```bash
   git add .
   git commit -m "Add completeCoboundaryTop optimization and enhanced timing"
   ```
   This saves the new code with all our improvements.

### Phase 2: Update `main` Branch with Timing (Keep Old Data Structures)

2. **Switch to `main` branch**:
   ```bash
   git checkout main
   ```

3. **Add timing enhancements to `main`** (but keep old data structures):
   - Update `readOFF()` and `readPLY()` with enhanced timing
   - Update Forman gradient timing
   - **Important**: Keep `adjRelations`, keep `buildDataStructure()` only (no parallel)
   - Add `USE_PARALLEL_BUILD` flag support (but it will always use sequential in old code)

4. **Create test program** (`test_forman_gradient.cpp`) - same as new code

5. **Update CMakeLists.txt** - add test program target

6. **Commit timing enhancements to `main`**:
   ```bash
   git add .
   git commit -m "Add enhanced timing instrumentation (old code with adjRelations)"
   ```

### Phase 3: Verify Both Branches

7. **Verify `main` branch**:
   - Has `adjRelations` in `simplicialcomplex.h`
   - Has enhanced timing output
   - Uses `buildDataStructure()` only
   - Test program compiles

8. **Verify `improve_iastar` branch**:
   ```bash
   git checkout improve_iastar
   ```
   - Has `completeCoboundaryTop` in `simplicialcomplex.h`
   - Has enhanced timing output
   - Supports both sequential and parallel builds
   - Test program compiles

### Phase 4: Compile Variants

9. **Compile old code variants from `main`**:
   ```bash
   git checkout main
   cd cpp/cmake-build-debug
   # Sequential (only option for old code)
   cmake ../.. -DCMAKE_BUILD_TYPE=Release -DUSE_PARALLEL_BUILD=OFF
   make xx_tts -j$(nproc)
   mv xx_tts xx_tts_seq_old
   ```

10. **Compile new code variants from `improve_iastar`**:
    ```bash
    git checkout improve_iastar
    cd cpp/cmake-build-debug
    # Sequential
    cmake ../.. -DCMAKE_BUILD_TYPE=Release -DUSE_PARALLEL_BUILD=OFF
    make xx_tts -j$(nproc)
    mv xx_tts xx_tts_seq_new
    
    # Parallel
    cmake ../.. -DCMAKE_BUILD_TYPE=Release -DUSE_PARALLEL_BUILD=ON
    make xx_tts -j$(nproc)
    mv xx_tts xx_tts_pa_new
    ```

## Detailed Implementation Steps

### Step 1: Commit Current Work

```bash
# Make sure we're on improve_iastar
git checkout improve_iastar

# Review what will be committed
git status

# Commit all changes
git add .
git commit -m "Add completeCoboundaryTop optimization and enhanced timing

- Remove adjRelations, use completeCoboundaryTop instead
- Add enhanced timing instrumentation
- Add test_forman_gradient program
- Add compile_variants.sh and run_tests.sh scripts
- Add comprehensive documentation"
```

### Step 2: Switch to Main and Add Timing

```bash
# Switch to main
git checkout main

# Verify we have old code (adjRelations)
grep -n "adjRelations" cpp/source/iastar/simplicialcomplex.h
```

### Step 3: Apply Timing Enhancements to Main

**Files to modify in `main`** (keeping old data structures):

1. **`cpp/source/iastar/io_functions.cpp`**:
   - Add enhanced timing for `readOFF()` and `readPLY()`
   - Keep using `buildDataStructure()` only (no parallel option)
   - Add timing output similar to new code

2. **`cpp/source/forman/formangradient.cpp`**:
   - Add enhanced timing for gradient encoding and filtration
   - Same timing structure as new code

3. **`cpp/source/projects/test_forman_gradient.cpp`**:
   - Copy from new code (same file)

4. **`cpp/CMakeLists.txt`**:
   - Add test program target (same as new code)

**Important**: Do NOT modify:
- `simplicialcomplex.h` - keep `adjRelations`
- `simplicialcomplex.cpp` - keep old `buildDataStructure()` logic
- Any data structure changes

### Step 4: Update Compilation Script

The `compile_variants.sh` script needs to be updated to:
- Use `main` branch for old code (not `old_adjRelations`)
- Use `improve_iastar` branch for new code
- Handle the fact that old code only supports sequential build

## Alternative Approach: Cherry-pick Timing Changes

Instead of manually applying timing changes, we could:

1. Create a separate branch from `main` with only timing changes
2. Cherry-pick timing-related commits (if we structure them separately)
3. This is more complex but keeps history cleaner

**Recommendation**: Manual application is simpler and clearer for this use case.

## Verification Checklist

After setup, verify:

- [ ] `main` branch:
  - [ ] Has `adjRelations` in `simplicialcomplex.h`
  - [ ] Does NOT have `completeCoboundaryTop`
  - [ ] Has enhanced timing output
  - [ ] Uses `buildDataStructure()` only
  - [ ] `test_forman_gradient` compiles

- [ ] `improve_iastar` branch:
  - [ ] Has `completeCoboundaryTop` in `simplicialcomplex.h`
  - [ ] Does NOT have `adjRelations`
  - [ ] Has enhanced timing output
  - [ ] Supports both sequential and parallel builds
  - [ ] `test_forman_gradient` compiles

- [ ] All 4 executables compile:
  - [ ] `xx_tts_seq_old` (from main, sequential)
  - [ ] `xx_tts_pa_old` (from main, sequential - same as seq_old)
  - [ ] `xx_tts_seq_new` (from improve_iastar, sequential)
  - [ ] `xx_tts_pa_new` (from improve_iastar, parallel)

## Notes

1. **Old code parallel variant**: Since old code only has `buildDataStructure()`, `xx_tts_pa_old` will be the same as `xx_tts_seq_old`. This is fine - we're mainly comparing old vs new.

2. **Timing consistency**: Both branches should have the same timing output format for easy comparison.

3. **Test program**: Same in both branches, so we can compare results directly.

## Next Steps After Setup

1. Compile all variants
2. Run tests on all test files
3. Compare performance (old vs new, sequential vs parallel)
4. Generate performance report

---

**Status**: Planning Phase
**Ready to Execute**: Yes, after review

