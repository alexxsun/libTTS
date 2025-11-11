# Timing Comparison Guide

## Quick Method: Use the Script

Run the automated comparison script:
```bash
cd /home/alex/Projects/libTTS_public
./compare_timing.sh
```

This will:
1. Save your current changes
2. Build and test the OLD version (with adjRelations)
3. Restore your changes
4. Build and test the NEW version (without adjRelations)
5. Show a comparison

## Manual Method (Step by Step)

If you prefer to do it manually:

### Step 1: Save Current Changes
```bash
cd /home/alex/Projects/libTTS_public
git stash push -m "New implementation" -- cpp/source/iastar/
```

### Step 2: Build Old Version
```bash
cd cpp/cmake-build-debug
cmake ../
make xx_tts -j4
```

### Step 3: Test Old Version
```bash
time ./xx_tts close_stems_3_a0.01.ply close_stems_3_locs.pts -tts 2>&1 | tee timing_old.log
```

**Save the output file:**
```bash
mv close_stems_3_a0.01_lbl.ply close_stems_3_a0.01_lbl_OLD.ply
```

### Step 4: Restore New Changes
```bash
cd /home/alex/Projects/libTTS_public
git stash pop
```

### Step 5: Build New Version
```bash
cd cpp/cmake-build-debug
cmake ../
make xx_tts -j4
```

### Step 6: Test New Version
```bash
time ./xx_tts close_stems_3_a0.01.ply close_stems_3_locs.pts -tts 2>&1 | tee timing_new.log
```

**Save the output file:**
```bash
mv close_stems_3_a0.01_lbl.ply close_stems_3_a0.01_lbl_NEW.ply
```

### Step 7: Compare Results

**Extract key timing metrics:**
```bash
echo "=== OLD VERSION ==="
grep -E "(build IA\*|Tops computed|grow time|label time)" timing_old.log

echo ""
echo "=== NEW VERSION ==="
grep -E "(build IA\*|Tops computed|grow time|label time)" timing_new.log

echo ""
echo "=== Complete mapping storage (NEW only) ==="
grep "Complete mapping storage" timing_new.log
```

**Compare overall time:**
```bash
echo "=== Overall Time Comparison ==="
echo "OLD:"
tail -3 timing_old.log | grep -E "real|user|sys"

echo ""
echo "NEW:"
tail -3 timing_new.log | grep -E "real|user|sys"
```

## What to Look For

### Key Metrics:

1. **Data Structure Building** (`ply build IA*`)
   - Expected: 20-40% faster in NEW version
   - Look for: Time reduction

2. **storeFullStar** (`Tops computed`)
   - Expected: Same or faster in NEW version
   - Look for: Time comparison

3. **Complete Mapping Storage** (NEW only)
   - Shows: Cost of storing complete mapping
   - Should be: Very small (< 0.1s typically)

4. **Overall Process Time**
   - Expected: Faster in NEW version
   - Look for: Total time saved

### Example Output Format:

**OLD VERSION:**
```
   ply build IA*: 3.579 s
Tops computed 5.123 s
grow time: 2.456 s
label time: 1.789 s

real    0m13.456s
user    0m12.234s
sys     0m1.123s
```

**NEW VERSION:**
```
      Complete mapping storage took: 0.123 seconds
   ply build IA*: 1.357 s
Tops computed 3.234 s
grow time: 2.456 s
label time: 1.789 s

real    0m8.234s
user    0m7.123s
sys     0m1.111s
```

**Improvement:**
- Build time: 3.579s → 1.357s (62% faster, saved 2.222s)
- storeFullStar: 5.123s → 3.234s (37% faster, saved 1.889s)
- Overall: 13.456s → 8.234s (39% faster, saved 5.222s)

## Troubleshooting

### If stash fails:
```bash
# Check what's in stash
git stash list

# If needed, commit current changes first
git add cpp/source/iastar/
git commit -m "New implementation without adjRelations"
```

### If build fails:
- Make sure you're in the correct directory
- Check that test files exist in build directory
- Verify CMake configuration

### If test files not found:
Copy them to the build directory:
```bash
cp /path/to/close_stems_3_a0.01.ply cpp/cmake-build-debug/
cp /path/to/close_stems_3_locs.pts cpp/cmake-build-debug/
```

## Alternative: Use Git Branches

If you prefer using branches:

```bash
# Create a branch for old version
git checkout -b old_version
git add cpp/source/iastar/
git commit -m "Old version with adjRelations"

# Build and test old version
cd cpp/cmake-build-debug
cmake ../ && make xx_tts -j4
time ./xx_tts close_stems_3_a0.01.ply close_stems_3_locs.pts -tts 2>&1 | tee timing_old.log

# Switch back to new version
cd /home/alex/Projects/libTTS_public
git checkout improve_iastar

# Build and test new version
cd cpp/cmake-build-debug
cmake ../ && make xx_tts -j4
time ./xx_tts close_stems_3_a0.01.ply close_stems_3_locs.pts -tts 2>&1 | tee timing_new.log
```

