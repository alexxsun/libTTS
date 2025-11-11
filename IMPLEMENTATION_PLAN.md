# Implementation Plan: Remove `adjRelations`

## Goal
Remove `adjRelations` data structure and replace it with `completeCoboundaryTop` while maintaining identical output.

## Test Setup
- **Test command**: `./xx_tts close_stems_3_a0.01.ply close_stems_3_locs.pts -tts`
- **Expected output**: `close_stems_3_a0.01_lbl.ply` (should match `old_close_stems_3_a0.01_lbl.ply`)
- **Build location**: `/home/alex/Projects/libTTS_public/cpp/cmake-build-debug/`

## Implementation Strategy: Incremental with Testing

### Phase 1: Add New Data Structure (Non-Breaking)
**Goal**: Add `completeCoboundaryTop` without removing anything yet.

#### Step 1.1: Add Data Structure Declaration
**File**: `cpp/source/iastar/simplicialcomplex.h`
- Add after line 29:
  ```cpp
  vector<map<int, set<int>>> completeCoboundaryTop;  // vertex → dimension → top simplex indices
  ```
- **Test**: Compile should succeed

#### Step 1.2: Initialize in Constructor
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- In `SimplicialComplex::SimplicialComplex()` (line 6), add:
  ```cpp
  completeCoboundaryTop.clear();
  ```
- **Test**: Compile and run test - should produce same output

#### Step 1.3: Store Complete Mapping in Sequential Build
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- In `buildDataStructure()`, after line 130 (after `incidentTop` is built):
  ```cpp
  // Store complete vertex-to-top-simplex mapping
  auto start_storage = std::chrono::high_resolution_clock::now();
  if (completeCoboundaryTop.size() <= vertices.size()) {
      completeCoboundaryTop.resize(vertices.size());
  }
  for (uint j = 0; j < vertices.size(); j++) {
      if (!incidentTop[j].empty()) {
          completeCoboundaryTop[j][dim].insert(incidentTop[j].begin(), incidentTop[j].end());
      }
  }
  auto end_storage = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> storage_elapsed = end_storage - start_storage;
  std::cout << "      Complete mapping storage took: " << storage_elapsed.count() << " seconds\n";
  ```
- **Test**: Compile and run test - should produce same output
- **Note**: Record timing for comparison

#### Step 1.4: Store Complete Mapping in Parallel Build
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- In `buildDataStructure_parallel()`, after line 299 (after `incidentTop` is built):
  ```cpp
  // Store complete vertex-to-top-simplex mapping
  auto start_storage = Clock::now();
  if (completeCoboundaryTop.size() <= vertices.size()) {
      completeCoboundaryTop.resize(vertices.size());
  }
  for (uint j = 0; j < vertices.size(); j++) {
      if (!incidentTop[j].empty()) {
          completeCoboundaryTop[j][dim].insert(incidentTop[j].begin(), incidentTop[j].end());
      }
  }
  auto end_storage = Clock::now();
  std::chrono::duration<double> storage_elapsed = end_storage - start_storage;
  std::cout << "      Complete mapping storage took: " << storage_elapsed.count() << " seconds\n";
  ```
- **Test**: Compile and run test - should produce same output
- **Note**: Record timing for comparison

**Phase 1 Verification**: 
- ✅ Code compiles
- ✅ Test produces identical output
- ✅ `completeCoboundaryTop` is populated

---

### Phase 2: Implement New `incidentCluster()` (Side-by-Side)
**Goal**: Create new version that uses `completeCoboundaryTop`, keep old one for comparison.

#### Step 2.1: Create New Function
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- Add new function `incidentCluster_new()` after line 478:
  ```cpp
  forward_list<explicitS> *SimplicialComplex::incidentCluster_new(explicitS vertex, explicitS topS) {
      forward_list<explicitS> *ret = new forward_list<explicitS>();
      set<explicitS> visited;
      queue<explicitS> adjacentSimplexes;
      
      adjacentSimplexes.push(topS);
      visited.insert(topS);
      
      while (!adjacentSimplexes.empty()) {
          explicitS current = adjacentSimplexes.front();
          TopSimplex &top = getTopSimplex(current);
          int dim = top.getDimension();
          
          // For each face (excluding faces containing the vertex)
          for (int i = 0; i < top.get_nVertices(); i++) {
              if (top.getVertexIndex(i) == vertex.getIndex()) {
                  continue; // Skip faces containing the vertex
              }
              
              // Build face vertices (all vertices except the i-th one)
              set<int> faceVertices;
              for (int j = 0; j < top.get_nVertices(); j++) {
                  if (j != i) {
                      faceVertices.insert(top.getVertexIndex(j));
                  }
              }
              
              // Find all top simplexes of same dimension incident to vertex that share this face
              auto &vertexMap = completeCoboundaryTop[vertex.getIndex()];
              auto dimIt = vertexMap.find(dim);
              if (dimIt == vertexMap.end()) continue;
              
              for (int candidateIdx : dimIt->second) {
                  explicitS candidate(dim, candidateIdx);
                  
                  if (visited.find(candidate) != visited.end()) {
                      continue; // Already processed
                  }
                  
                  TopSimplex &candidateTop = getTopSimplex(candidate);
                  vector<int> &candidateVertices = candidateTop.getVertices();
                  
                  // Check if candidate shares the face (has all face vertices)
                  bool sharesFace = true;
                  for (int fv : faceVertices) {
                      bool found = false;
                      for (int cv : candidateVertices) {
                          if (cv == fv) {
                              found = true;
                              break;
                          }
                      }
                      if (!found) {
                          sharesFace = false;
                          break;
                      }
                  }
                  
                  if (sharesFace) {
                      visited.insert(candidate);
                      adjacentSimplexes.push(candidate);
                  }
              }
          }
          
          adjacentSimplexes.pop();
      }
      
      ret->insert_after(ret->before_begin(), visited.begin(), visited.end());
      return ret;
  }
  ```

#### Step 2.2: Add Declaration
**File**: `cpp/source/iastar/simplicialcomplex.h`
- Add in protected section (around line 207):
  ```cpp
  forward_list<explicitS> *incidentCluster_new(explicitS vertex, explicitS topS);
  ```

#### Step 2.3: Test New Function (Temporary Switch)
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- In `topStar()` function, temporarily replace line 502:
  ```cpp
  // OLD: forward_list<explicitS> *adjs = incidentCluster(vertex, simplexes[i]);
  forward_list<explicitS> *adjs = incidentCluster_new(vertex, simplexes[i]);
  ```
- Also in `topStar(vertex, dimension)` at line 524:
  ```cpp
  // OLD: forward_list<explicitS> *adjs = incidentCluster(vertex, simplexes[i]);
  forward_list<explicitS> *adjs = incidentCluster_new(vertex, simplexes[i]);
  ```
- **Test**: Compile and run test - should produce **identical** output
- **If output differs**: Debug and fix `incidentCluster_new()` until output matches

**Phase 2 Verification**:
- ✅ New function compiles
- ✅ Test produces **identical** output to old version
- ✅ Both functions produce same results

---

### Phase 3: Switch to New Implementation
**Goal**: Replace old `incidentCluster()` with new one.

#### Step 3.1: Replace Function Implementation
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- Replace entire `incidentCluster()` function (lines 445-478) with `incidentCluster_new()` code
- Remove `_new` suffix from function name

#### Step 3.2: Update Function Calls
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- In `buildDataStructure()` line 135, ensure it calls `incidentCluster()` (should already)
- In `buildDataStructure_parallel()` line 328, ensure it calls `incidentCluster()` (should already)

#### Step 3.3: Remove Temporary Declaration
**File**: `cpp/source/iastar/simplicialcomplex.h`
- Remove `incidentCluster_new()` declaration

**Phase 3 Verification**:
- ✅ Code compiles
- ✅ Test produces identical output
- ✅ Old `incidentCluster()` is replaced

---

### Phase 4: Remove `adjRelations` Building Code
**Goal**: Remove expensive adjacency building while keeping functionality.

#### Step 4.1: Remove from Sequential Build
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- Remove lines 41-114 (adjacency relations building phase)
- Keep the vertex incidence building (lines 119-145)

#### Step 4.2: Remove from Parallel Build
**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- Remove lines 176-284 (adjacency relations building phase)
- Keep the vertex incidence building (lines 290-350)

#### Step 4.3: Test
- **Test**: Compile and run test - should produce identical output

**Phase 4 Verification**:
- ✅ Code compiles
- ✅ Test produces identical output
- ✅ Build time should be faster

---

### Phase 5: Remove `adjRelations` Data Structure
**Goal**: Clean up unused code.

#### Step 5.1: Remove Declaration
**File**: `cpp/source/iastar/simplicialcomplex.h`
- Remove line 29: `vector<forward_list<int> > adjRelations;`

#### Step 5.2: Remove `topAdjacent()` Function
**File**: `cpp/source/iastar/simplicialcomplex.h`
- Remove declaration at line 120

**File**: `cpp/source/iastar/simplicialcomplex.cpp`
- Remove function implementation (lines 428-443)

#### Step 5.3: Remove from I/O
**File**: `cpp/source/iastar/io_functions.cpp`
- Remove `adjRelations` loading code in `readIA()` (lines 701-713)
- **Note**: If `saveIA()` exists and saves `adjRelations`, we have two options:
  1. Remove saving code (if files are regenerated)
  2. Keep compatibility by saving empty `adjRelations` (if old files need to be readable)
- **Check**: Search for `saveIA` implementation to see if it saves `adjRelations`

#### Step 5.4: Final Test
- **Test**: Compile and run test - should produce identical output

**Phase 5 Verification**:
- ✅ Code compiles
- ✅ Test produces identical output
- ✅ No references to `adjRelations` remain

---

## Testing Checklist

After each phase, verify:

1. **Compilation**: `cd /home/alex/Projects/libTTS_public/cpp/cmake-build-debug && cmake ../ && make xx_tts -j4`
2. **Execution with Timing**: `./xx_tts close_stems_3_a0.01.ply close_stems_3_locs.pts -tts`
3. **Output Comparison**: Compare `close_stems_3_a0.01_lbl.ply` with `old_close_stems_3_a0.01_lbl.ply`
   ```bash
   # Quick check: file size and basic properties
   ls -lh close_stems_3_a0.01_lbl.ply old_close_stems_3_a0.01_lbl.ply
   
   # Detailed comparison (if you have a comparison tool)
   # diff or custom comparison script
   ```
4. **Performance Timing**: Record and compare:
   - Data structure building time (from console output)
   - Overall process time (from console output)
   - See "Timing Measurements" section below

## Timing Measurements

### Current Timing Points

The code already has some timing in place:
- `buildDataStructure()` timing in `readPLY()` (line 203-208 in `io_functions.cpp`)
- `buildDataStructure_parallel()` has `std::chrono` timing (lines 152-359)
- `storeFullStar()` timing in `computeFormanGradient()` (line 110-116)

### Adding Detailed Timing

We'll add timing to measure:
1. **Adjacency building phase** (before removal) - to see what we're eliminating
2. **Complete mapping storage** (new) - to see the cost
3. **Overall build time** - to see total improvement

### Step-by-Step Timing Additions

#### Phase 1: Add Timing to Current Code (Baseline)
**File**: `cpp/source/iastar/simplicialcomplex.cpp`

**In `buildDataStructure()`** (sequential version):
- Add timing around adjacency building (lines 41-114):
  ```cpp
  #include <chrono>
  using Clock = std::chrono::high_resolution_clock;
  
  // Before line 41:
  auto start_adj = Clock::now();
  
  // ... adjacency building code (lines 41-114) ...
  
  // After line 114:
  auto end_adj = Clock::now();
  std::chrono::duration<double> adj_elapsed = end_adj - start_adj;
  std::cout << "      Adjacency building took: " << adj_elapsed.count() << " seconds\n";
  ```

**In `buildDataStructure_parallel()`**:
- The timing is already there (lines 176, 286-288), but we'll enhance it:
  ```cpp
  // Already exists at line 176: auto st_adj_phase = Clock::now();
  // Already exists at line 286-288: timing output
  // We'll keep this and add comparison
  ```

#### Phase 2-5: Track Timing Improvements

After each phase, the console output will show:
- **Before**: Adjacency building time (from Phase 1)
- **After Phase 4**: No adjacency building (removed), but complete mapping storage time
- **Comparison**: Time saved

### Expected Timing Output Format

**Before (with adjRelations)**:
```
dim: 2, simplex #: 12345, faces #: 37035
      Adjacency building took: 2.345 seconds
      Phase 0, adj phase took: 2.345 seconds
      Phase 0, vertex addPartialCoboundaryTop took: 1.234 seconds
Phase 0 total took: 3.579 seconds

   ply build IA*: 3.579 s
Tops computed 5.123 s
grow time: 2.456 s
label time: 1.789 s
```

**After (without adjRelations)**:
```
dim: 2, simplex #: 12345, faces #: 37035
      Complete mapping storage took: 0.123 seconds
      Phase 0, vertex addPartialCoboundaryTop took: 1.234 seconds
Phase 0 total took: 1.357 seconds

   ply build IA*: 1.357 s  [IMPROVEMENT: 2.222s saved, 62% faster]
Tops computed 3.234 s  [IMPROVEMENT: 1.889s saved, 37% faster]
grow time: 2.456 s  [Same]
label time: 1.789 s  [Same]
```

### Timing Measurement Script

Create a simple script to extract and compare timings:

```bash
#!/bin/bash
# compare_timing.sh

echo "=== Timing Comparison ==="
echo ""
echo "BEFORE (with adjRelations):"
grep -E "(build IA\*|Tops computed|grow time|label time)" old_output.log

echo ""
echo "AFTER (without adjRelations):"
grep -E "(build IA\*|Tops computed|grow time|label time)" new_output.log

echo ""
echo "=== Improvement Summary ==="
# Calculate improvements (manual or with script)
```

### Manual Timing Collection

Run the test and collect output:
```bash
cd /home/alex/Projects/libTTS_public/cpp/cmake-build-debug/

# Capture full output with timing
time ./xx_tts close_stems_3_a0.01.ply close_stems_3_locs.pts -tts 2>&1 | tee timing_output.log
```

The `time` command will show:
- **Real time**: Wall-clock time (what you experience)
- **User time**: CPU time in user mode
- **Sys time**: CPU time in system mode

Key metrics to record from console output:
1. **`ply build IA*`**: Data structure building time
2. **`Tops computed`**: `storeFullStar()` time (includes `incidentCluster()` calls)
3. **`grow time`**: Tree growing time (should be same)
4. **`label time`**: Labeling time (should be same)
5. **Total time** (from `time` command): Overall process time

### Optional: Add Overall Timing to Code

If you want more precise overall timing, we can add it to `extract_single_trees()`:

**File**: `cpp/source/projects/xx_tts.cpp`
- Add at the start of `extract_single_trees()` (line 93):
  ```cpp
  IO_Timer total_timer;
  total_timer.start();
  ```
- Add at the end (before return):
  ```cpp
  total_timer.stop();
  cout << "\n=== TOTAL PROCESS TIME: " << total_timer.getElapsedTimeInSec() << " s ===\n";
  ```

### Performance Targets

Based on analysis:
- **Data structure building**: 20-40% faster (eliminate sorting)
- **`storeFullStar()`**: Potentially faster (faster `incidentCluster()`)
- **Overall process**: Should see improvement proportional to build time savings

## Rollback Strategy

If any phase fails:
1. **Git commit** before starting: `git commit -am "Before adjRelations removal"`
2. **After each phase**: `git commit -am "Phase X complete"`
3. **If issues**: `git reset --hard HEAD` to rollback to last good commit

## Expected Outcomes

### Performance Improvements
- **Build time**: 20-40% faster (no sorting phase)
- **Memory**: Slightly different (eliminate `adjRelations`, add `completeCoboundaryTop`)
- **Runtime**: Same (cache still used)

### Code Changes
- **Lines removed**: ~230
- **Lines added**: ~50
- **Net**: ~180 lines removed

### Timing Results Template

After implementation, document results:

```
=== Performance Results ===

Test file: close_stems_3_a0.01.ply

BEFORE (with adjRelations):
- Data structure building: X.XXX seconds
- storeFullStar(): X.XXX seconds
- Overall process: X.XXX seconds

AFTER (without adjRelations):
- Data structure building: X.XXX seconds (Y% faster)
- storeFullStar(): X.XXX seconds (Y% faster)
- Overall process: X.XXX seconds (Y% faster)

Improvements:
- Build time saved: X.XXX seconds
- Total time saved: X.XXX seconds
- Speedup: Y%
```

## Risk Assessment

### Low Risk
- Phase 1: Adding new data structure (non-breaking)
- Phase 2: Side-by-side implementation (can compare)

### Medium Risk
- Phase 3: Replacing function (but already tested)
- Phase 4: Removing building code (but functionality preserved)

### Higher Risk
- Phase 5: Removing data structure (final cleanup, but should be safe)

## Notes

1. **Parallel Build**: Make sure both sequential and parallel builds work
2. **I/O Compatibility**: If files are saved/loaded, ensure compatibility (or version bump)
3. **Edge Cases**: Test with different mesh sizes and topologies
4. **Debugging**: Add temporary debug output if needed to compare intermediate results


