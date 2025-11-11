# Test Results Interpretation Guide

## Issue: Empty Log Files

The log files in your test results are empty (0 bytes). This is because:

1. **`xx_tts` requires specific command-line arguments** - it doesn't work with just a filename
2. **For `.off` files**: You should use `test_forman_gradient` instead of `xx_tts`
3. **For `.ply` files**: `xx_tts` needs the `-tts` argument with a trunk file

## Solution: Run Tests Manually

### For `.off` Files (Forman Gradient Testing)

Use `test_forman_gradient` program:

```bash
cd cpp/cmake-build-debug

# For each variant, compile test_forman_gradient on the appropriate branch
# Then run:

./test_forman_gradient close_stems_3_a0.010.off 3 > close_stems_3_seq_old.log 2>&1
./test_forman_gradient close_stems_3_a0.010.off 3 > close_stems_3_seq_new.log 2>&1
./test_forman_gradient close_stems_3_a0.010.off 3 > close_stems_3_pa_new.log 2>&1
# (pa_old is same as seq_old for old code)
```

### For `.ply` Files (Complete Workflow)

Use `xx_tts` with proper arguments:

```bash
cd cpp/cmake-build-debug

# Run with -tts argument (requires trunk file)
./xx_tts_seq_old close_stems_3_a0.01.ply close_stems_3_locs.pts -tts > ply_seq_old.log 2>&1
./xx_tts_seq_new close_stems_3_a0.01.ply close_stems_3_locs.pts -tts > ply_seq_new.log 2>&1
./xx_tts_pa_new close_stems_3_a0.01.ply close_stems_3_locs.pts -tts > ply_pa_new.log 2>&1
```

## Updated Script

I've updated `run_tests.sh` to:
1. ✅ Use `tee` for output redirection (as you suggested)
2. ✅ Use proper arguments for `.ply` files (`-tts` with trunk file)
3. ✅ Better handle `test_forman_gradient` for `.off` files
4. ✅ Check if output files are empty and warn

## Next Steps

1. **Re-run tests with updated script**:
   ```bash
   ./run_tests.sh cpp/cmake-build-debug
   ```

2. **Or run tests manually** (recommended for now):
   ```bash
   # Create a directory for results
   mkdir -p cpp/cmake-build-debug/manual_test_results
   cd cpp/cmake-build-debug
   
   # Test .off files with test_forman_gradient
   ./test_forman_gradient close_stems_3_a0.010.off 3 | tee ../manual_test_results/close_stems_3_seq_new.log
   
   # Test .ply files with xx_tts
   ./xx_tts_seq_new close_stems_3_a0.01.ply close_stems_3_locs.pts -tts | tee ../manual_test_results/ply_seq_new.log
   ```

3. **Analyze results**:
   ```bash
   python3 interpret_results.py cpp/cmake-build-debug/manual_test_results
   ```

## What to Extract from Logs

Once you have log files with content, look for:

1. **File Statistics**:
   - `Vertices: X`
   - `Total top simplexes: X`
   - `Complex dimension: X`

2. **Timing Metrics**:
   - `ply/off read file I/O: X.XXX s`
   - `ply/off read vertices: X.XXX s`
   - `ply/off read cells: X.XXX s`
   - `ply/off build IA* (sequential/parallel): X.XXX s`
   - `ply/off total reading time: X.XXX s`
   - `Gradient encoding time: X.XXX s`
   - `Filtration computation time: X.XXX s`
   - `Total Forman gradient time: X.XXX s`
   - `Total time: X.XXX s`

3. **Calculate Speedups**:
   - Build IA* speedup = `time_old / time_new`
   - Forman gradient speedup = `time_old / time_new`
   - Total speedup = `time_old / time_new`

## Quick Analysis Command

Once you have log files with content:

```bash
# Analyze any test results directory
python3 interpret_results.py <test_results_dir>

# Or use the helper script
./analyze_test_results.sh cpp/cmake-build-debug
```

---

**Status**: Script updated with `tee` and proper argument handling. Ready to re-run tests.

