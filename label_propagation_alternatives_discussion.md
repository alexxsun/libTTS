# Label Propagation Methods: Alternatives Discussion

## ⚠️ PRIMARY CHALLENGE: Very Dense Point Clouds

**The biggest challenge is the extremely high point density:**
- Input point clouds contain **millions of points**
- **Time complexity is very high** - this is the primary bottleneck
- Current region growing method: O(M × N × log N) where M = labeled points, N = total points
- With millions of points, this becomes prohibitively slow
- **Performance is the critical constraint**, not just accuracy

**This document prioritizes methods that scale efficiently to dense point clouds.**

---

## Problem Context

**Current Situation:**
- You have **already labeled points**: tree trunks and main branches
- You need to label **remaining points**: leaves and tiny branches
- **Multiple trees** in the scene (need to maintain tree identity)
- Goal: Create more complete tree point clouds

**Key Challenges:**

1. **⚠️ VERY DENSE POINT CLOUDS - PRIMARY CHALLENGE**: 
   - Input points are **extremely dense** (millions of points)
   - **Time complexity is very high** - this is the biggest bottleneck
   - Need methods that scale efficiently to large point clouds
   - Memory usage must be manageable
   - KDTree queries become expensive with millions of points

2. **Multiple trees**: Need to ensure leaves/branches are labeled with correct tree ID

3. **Complex geometry**: Leaves and tiny branches may be disconnected or far from labeled points

4. **Density variations**: Leaves may be sparse, branches may be dense

5. **Occlusion**: Some parts may be hidden or disconnected

---

## Original Method: Z-Ordered Region Growing

**How it works:**
- Starts from labeled seed points (trunks/branches)
- Grows labels outward within a search radius
- Processes points in Z-order (bottom-up)
- Maintains tree identity through propagation

**Limitations:**
- ⚠️ **KDTree queries are expensive** with millions of points (O(N log N) per query)
- ⚠️ **Sequential processing** - each point queries neighbors individually
- ⚠️ **Memory intensive** - needs to store all points and labels
- May not reach disconnected components (leaves far from branches)
- Fixed search radius may miss sparse regions
- Sequential processing can be slow with many trees

**Time Complexity:**
- For N points, M labeled seeds: O(M × N × log N) in worst case
- With millions of points, this becomes prohibitively slow

---

## Implemented Methods

### Method 1: Z-Layer Batched Region Growing

**Concept**: Process points in Z-layers/batches while maintaining global Z-ordering.

**Status**: ✅ Implemented

**Code Location**: `libtts/label_propagation.py::label_points_region_growing_layered()`

**Performance**: 3-5x speedup through batch processing

**When to use**: When Z-ordering is critical and you have moderate point densities

**Limitations**: Still slower than distance-based methods

---

### Method 2: Distance-Based Label Propagation ⭐

**Concept**: For each unlabeled point, find the closest labeled point within a maximum distance and assign that label. One-shot operation (no iterative expansion).

**Status**: ✅ Implemented

**Code Location**: `libtts/label_propagation.py::label_points_distance_based()`

**Key Innovation**:
- Builds KDTree only on the much smaller set of *labeled* points (not all points)
- Single vectorized query for all unlabeled points
- No iterative expansion needed
- Fully parallelizable

**Time Complexity**: O(N × log M) where N = unlabeled points, M = labeled points (M << N)

**Performance**:
- ⚠️ **Issue Found**: Initially only labeled seed points (no propagation)
- ✅ **Fixed**: Added seed label projection step
- **Current**: ~520s for 10M points (but only 9.7% coverage - needs investigation)

**Advantages**:
- ✅ Fastest method for dense clouds
- ✅ Scales well (depends on labeled points, not total)
- ✅ Fully parallelizable
- ✅ Simple to implement

**Limitations**:
- ⚠️ Coverage issue (only 9.7% - needs investigation)
- May assign wrong tree if trees are close
- Doesn't consider connectivity

**Usage**:
```python
run_label_propagation(
    infile="unlabeled_cloud.ply",
    labeled_file="seed_points.ply",
    method='distance_based',
    max_distance=0.25,
    n_jobs=4,
    out_file="output.ply"
)
```

---

### Method 3: Iterative Distance-Based Wave Propagation ⭐⭐⭐ **RECOMMENDED**

**Concept**: Propagate labels in small iterative "waves", where each wave expands labels by a small distance (`wave_distance`) from currently labeled points. This creates gradual expansion similar to region growing but using distance queries.

**Status**: ✅ Implemented

**Code Location**: `libtts/label_propagation.py::label_points_iterative_distance_based()`

**Algorithm Overview**:

**Each iteration (wave):**
1. Build KDTree on currently labeled points
2. Query all unlabeled points to find labeled neighbors within `wave_distance` (e.g., 0.05m)
3. Apply selected strategy to handle multiple labels (if multiple labeled points found)
4. Update labels for newly labeled points
5. Check stopping criteria

**Key Parameters:**
- `wave_distance`: Small distance per wave (default: 0.05m)
- `max_iterations`: Safety limit (default: 5)
- `min_new_points`: Stop if fewer than this many points labeled in a wave (default: 10)
- `multiple_label_strategy`: How to handle multiple labeled points within threshold

**Stopping Criteria**:
1. No new points labeled (most common)
2. Maximum iterations reached
3. Minimum new points threshold
4. All points labeled (early exit)

**Multiple Label Strategies**:

| Strategy | Description | Performance | Status |
|----------|-------------|-------------|--------|
| **'closest'** ⭐ | Use closest labeled point | ~23s for 10M points | ✅ **Recommended** |
| **'majority'** | Use most common label | ~225s for 10M points | ⚠️ Too slow |
| **'weighted'** | Weight by inverse distance | N/A | ⚠️ Not optimized |
| **'hybrid'** | Majority → Weighted → Closest | ~1500s for 10M points | ⚠️ Very slow |
| **'same_only'** | Only if all neighbors same | N/A | ⚠️ Not optimized |

**Performance**:
- **'closest' strategy**: ~23s for 10M points, 97.1% coverage
- **'majority' strategy**: ~225s for 10M points, 97.1% coverage
- **Time Complexity**: O(W × N × log M) where W = number of waves, N = unlabeled points, M = labeled points

**Advantages**:
- ✅ Gradual propagation (similar to region growing)
- ✅ Small buffer per step (0.05m per wave)
- ✅ Handles multiple labels (multiple strategies)
- ✅ Clear stopping criteria
- ✅ Fast (vectorized KDTree queries)
- ✅ Maintains tree identity

**Usage**:
```python
# Recommended: Use 'closest' strategy (default)
run_label_propagation(
    infile="unlabeled_cloud.ply",
    labeled_file="seed_points.ply",
    method='iterative_distance_based',
    wave_distance=0.05,
    max_iterations=5,
    multiple_label_strategy='closest',  # Default and recommended
    out_file="output.ply"
)
```

**See**: `iterative_label_propagation_performance.md` for detailed performance analysis

---

## Alternative Methods (Not Implemented)

### Method 4: Multi-Scale Region Growing

**Concept**: Use different search radii at different scales to capture both nearby and distant points.

**Performance**: ⚠️ Still expensive - multiple passes make it slower than single-pass

**Status**: ❌ Not implemented - not recommended for dense clouds

---

### Method 5: Graph-Based Propagation

**Concept**: Build a connectivity graph and propagate labels through graph edges.

**Performance**: ❌ Very expensive - O(N²) operations, not scalable

**Status**: ❌ Not implemented - not recommended for dense clouds

---

### Method 6: Hierarchical Clustering + Label Assignment

**Concept**: Cluster points hierarchically, then assign labels to clusters based on labeled points within them.

**Performance**: ❌ Very slow - clustering millions of points is expensive

**Status**: ❌ Not implemented - not recommended for dense clouds

---

### Method 7: Coarse-to-Fine with Spatial Downsampling

**Concept**: Downsample point cloud, label at coarse scale, then refine at fine scale.

**Performance**: ✅ Very fast - work with 10% of points first

**Status**: ❌ Not implemented - could be useful for extremely dense clouds (>10M points)

**When to use**: When you have extremely dense clouds (>10M points) and memory is limited

---

### Method 8: Hybrid: Region Growing + Distance Fallback

**Concept**: Combine region growing (for connectivity) with distance-based (for disconnected components).

**Performance**: ⚠️ Medium - first phase is slow, second phase is fast

**Status**: ❌ Not implemented - could be useful for balanced approach

**When to use**: When you need both connectivity preservation and disconnected component handling

---

## Comparison Table

| Method | Status | Speed (10M pts) | Coverage | Scalability | Best For |
|--------|--------|-----------------|----------|-------------|----------|
| **Z-Ordered Region Growing** | ✅ Original | Slow | Good | Poor | Connected components |
| **Z-Layer Batched** | ✅ Implemented | Moderate | Good | Medium | Z-ordering critical |
| **Distance-Based** | ✅ Implemented | ~520s | 9.7% ⚠️ | Excellent | Dense clouds (needs fix) |
| **Iterative Distance-Based** | ✅ Implemented | **~23s** | **97.1%** | Excellent | **Dense clouds** ⭐⭐⭐ |
| **Multi-Scale** | ❌ Not implemented | Slow | Good | Poor | Multiple scales |
| **Graph-Based** | ❌ Not implemented | Very slow | Good | Very poor | Not for dense clouds |
| **Clustering** | ❌ Not implemented | Very slow | Good | Very poor | Not for dense clouds |
| **Coarse-to-Fine** | ❌ Not implemented | Very fast | Good | Excellent | Very dense clouds |
| **Hybrid** | ❌ Not implemented | Medium | Good | Medium | Balanced approach |

---

## Recommendations

### ⭐⭐⭐ **Primary Recommendation: Iterative Distance-Based with 'closest' Strategy**

**Why**:
- ✅ **Fastest**: ~23s for 10M points (5 waves)
- ✅ **Excellent coverage**: 97.1% of points labeled
- ✅ **Fully optimized**: Vectorized, no Python loops
- ✅ **Recommended as default**: Use for all new code

**Usage**:
```python
# Default - uses 'closest' strategy
run_label_propagation(
    infile="unlabeled_cloud.ply",
    labeled_file="seed_points.ply",
    method='iterative_distance_based',
    out_file="output.ply"
)
```

### ⚠️ **Distance-Based Method: Needs Investigation**

**Current Status**:
- Only labels 9.7% of points (coverage issue)
- Needs investigation to fix propagation

**Future**: Once fixed, could be very fast alternative

### ❌ **Not Recommended for Dense Clouds**

- Graph-Based: Too slow (O(N²) complexity)
- Hierarchical Clustering: Too slow (O(N log N) to O(N²))
- Pure Region Growing: Too slow for millions of points
- Multi-Scale Region Growing: Multiple passes make it even slower

---

## Implementation Status

### ✅ Implemented and Recommended

1. **Iterative Distance-Based Wave Propagation** (`iterative_distance_based`)
   - Default method
   - Use 'closest' strategy
   - Excellent performance and coverage

### ✅ Implemented but Needs Work

2. **Distance-Based Label Propagation** (`distance_based`)
   - Coverage issue (only 9.7%)
   - Needs investigation

3. **Z-Layer Batched Region Growing** (`region_growing_layered`)
   - Works but slower than iterative method
   - Use when Z-ordering is critical

### ❌ Not Implemented (Future Options)

4. **Coarse-to-Fine with Downsampling**
   - Could be useful for extremely dense clouds (>10M points)
   - Low priority - iterative method already handles dense clouds well

5. **Hybrid Method**
   - Could be useful for balanced approach
   - Low priority - iterative method already provides good results

---

## Key Insights

### The Bottleneck

**The bottleneck is KDTree queries on millions of points.**

**Solution**: Build KDTree on labeled points only (much smaller!)
- **Result**: O(N × log M) instead of O(N × log N) where M << N
- **Example**: 1M points, 10K labeled → 100x fewer points in KDTree

### Why Iterative Method Works Well

1. **Gradual expansion**: Small waves (0.05m) prevent over-labeling
2. **Vectorized queries**: Fast KDTree queries on labeled points
3. **Vectorized strategy**: 'closest' strategy fully vectorized (663x speedup)
4. **Stopping criteria**: Clear conditions prevent infinite loops
5. **Excellent coverage**: 97.1% of points labeled

### Performance Comparison

| Method | Wave 1 Time | Total Time | Coverage | Recommendation |
|--------|-------------|------------|----------|----------------|
| **Iterative (closest)** | 9.51s | **~23s** | **97.1%** | ✅ **Use this** |
| **Iterative (majority)** | 205.67s | ~225s | 97.1% | ⚠️ Too slow |
| **Distance-Based** | N/A | ~520s | 9.7% | ⚠️ Needs fix |
| **Region Growing** | N/A | Very slow | Good | ❌ Not recommended |

---

## Future Work

### High Priority

1. **Fix Distance-Based Method Coverage Issue**
   - Currently only labels 9.7% of points
   - Investigate why propagation isn't working
   - Once fixed, could be fast alternative

### Medium Priority

1. **Optimize 'majority' Strategy Further** (if needed)
   - Currently 10x slower than 'closest'
   - May not be worth it - 'closest' is already excellent

2. **Vectorize 'hybrid' and 'weighted' Strategies** (if needed)
   - Currently very slow
   - Low priority - 'closest' may be sufficient

### Low Priority

1. **Implement Coarse-to-Fine Method**
   - For extremely dense clouds (>10M points)
   - Low priority - iterative method already handles dense clouds well

2. **Implement Hybrid Method**
   - For balanced approach
   - Low priority - iterative method already provides good results

---

## Conclusion

The **Iterative Distance-Based Wave Propagation** method with **'closest' strategy** is the clear winner:
- ✅ **Fastest**: ~23s for 10M points (5 waves)
- ✅ **Excellent coverage**: 97.1% of points labeled
- ✅ **Fully optimized**: Vectorized, no Python loops
- ✅ **Recommended as default**: Use for all new code

Other methods are available but not recommended until further optimization or investigation.

---

*Last Updated: 2024*
*Recommended Method: `iterative_distance_based` with `multiple_label_strategy='closest'`*
*See `iterative_label_propagation_performance.md` for detailed performance analysis*
