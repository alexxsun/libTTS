# Optimization Plan: Simplified IA* with Complete Vertex-to-Top Mapping

## Current State

✅ **Already Implemented:**
- `completeCoboundaryTop` stores ALL top simplexes per vertex per dimension
- No `adjRelations` building (removed)
- OpenMP parallelization in place

## Optimization Opportunities

### 1. ✅ **Optimize `topStar()` - Direct Lookup (HIGH IMPACT)**

**Current Implementation:**
- Uses `partialCoboundaryTop` (only representatives)
- Calls `incidentCluster()` multiple times to expand clusters
- Slow: O(k * m * d) where k = clusters, m = candidates, d = dimension

**Optimized Implementation:**
- Directly use `completeCoboundaryTop[vertex]` 
- No need for `incidentCluster()` calls
- Fast: O(n) where n = number of top simplexes (already stored)

**Speedup Expected:** 10-100x faster for `topStar()` calls

### 2. ✅ **Simplify `storeFullStar()` - Direct Copy**

**Current Implementation:**
- Calls `topStar()` which calls `incidentCluster()` multiple times
- Expensive cluster expansion

**Optimized Implementation:**
- Directly copy from `completeCoboundaryTop`
- No cluster expansion needed
- Just flatten the map into vectors

**Speedup Expected:** 10-50x faster for `storeFullStar()`

### 3. ⚠️ **Optimize `incidentCluster()` - Still Needed for Connectivity**

**Current Implementation:**
- Uses BFS traversal to find connected components
- Still needed for finding clusters of adjacent top simplexes

**Note:** We can't eliminate `incidentCluster()` entirely because it finds **connected components** (adjacent top simplexes sharing faces), not just all incident top simplexes. However, we can optimize it using `completeCoboundaryTop`.

## Implementation Strategy

### Phase 1: Optimize `topStar()` (Simple & High Impact)

Replace the current `topStar()` implementation to directly use `completeCoboundaryTop`:

```cpp
vector<explicitS> *SimplicialComplex::topStar(const explicitS &vertex) {
    assert(vertex.getDim() == 0);
    
    // Check cache first
    if (topPerVertex.size() > vertex.getIndex() && topPerVertex[vertex.getIndex()].size() != 0) {
        return new vector<explicitS>(topPerVertex[vertex.getIndex()]);
    }
    
    // Direct lookup from completeCoboundaryTop - MUCH FASTER!
    vector<explicitS> ret;
    int vertexIdx = vertex.getIndex();
    
    if (vertexIdx < completeCoboundaryTop.size()) {
        auto &vertexMap = completeCoboundaryTop[vertexIdx];
        for (const auto &dimPair : vertexMap) {
            int dim = dimPair.first;
            for (int topIdx : dimPair.second) {
                ret.push_back(explicitS(dim, topIdx));
            }
        }
    }
    
    return new vector<explicitS>(ret);
}
```

### Phase 2: Optimize `storeFullStar()` (Very Simple)

Since `topStar()` is now fast, `storeFullStar()` will automatically be faster. But we can make it even faster by directly copying:

```cpp
void SimplicialComplex::storeFullStar() {
    topPerVertex = vector<vector<explicitS> >(getVerticesNum());
    
#pragma omp parallel for schedule(dynamic) num_threads(6)
    for (int i = 0; i < getVerticesNum(); i++) {
        vector<explicitS> tops;
        
        if (i < completeCoboundaryTop.size()) {
            auto &vertexMap = completeCoboundaryTop[i];
            for (const auto &dimPair : vertexMap) {
                int dim = dimPair.first;
                for (int topIdx : dimPair.second) {
                    tops.push_back(explicitS(dim, topIdx));
                }
            }
        }
        
        sort(tops.begin(), tops.end());
        topPerVertex[i] = tops;
    }
}
```

### Phase 3: Keep `incidentCluster()` Optimized (Already Done)

The current `incidentCluster()` already uses `completeCoboundaryTop` efficiently. It's still needed for finding connected components.

## Expected Performance Improvements

### Before Optimization:
- `topStar()`: ~0.001-0.01s per call (with cluster expansion)
- `storeFullStar()`: ~5s (sequential cluster expansion for all vertices)

### After Optimization:
- `topStar()`: ~0.0001s per call (direct lookup)
- `storeFullStar()`: ~0.1-0.5s (direct copy, parallelized)

**Overall Speedup:**
- `topStar()`: **10-100x faster**
- `storeFullStar()`: **10-50x faster**
- Total initialization: **2-5x faster overall**

## Memory Trade-off

**Memory Usage:**
- `completeCoboundaryTop`: Stores all top simplexes per vertex
- Memory overhead: ~O(V * avg_tops_per_vertex * sizeof(int))
- For typical meshes: 10-50MB additional memory
- **Worth it** for 10-100x speedup!

## Implementation Priority

1. ✅ **HIGH**: Optimize `topStar()` - simple change, huge impact
2. ✅ **HIGH**: Optimize `storeFullStar()` - simple change, huge impact  
3. ⚠️ **MEDIUM**: Further optimize `incidentCluster()` if needed (already optimized)

## Code Changes Summary

**Files to Modify:**
- `cpp/source/iastar/simplicialcomplex.cpp`:
  - `topStar(explicitS)` - replace with direct lookup
  - `topStar(explicitS, int)` - replace with direct lookup
  - `storeFullStar()` - optimize to direct copy

**Lines Changed:** ~50 lines
**Complexity:** Low (straightforward refactoring)
**Risk:** Low (well-defined behavior)

