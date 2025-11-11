#!/usr/bin/env python3
"""
Interpret and analyze test results from test_results directory.

Usage:
    python interpret_results.py <test_results_dir>
    
Example:
    python interpret_results.py cpp/cmake-build-debug/test_results_20251111_114634
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict

def parse_timing_log(log_file):
    """Parse timing information from a log file."""
    data = {}
    
    if not os.path.exists(log_file) or os.path.getsize(log_file) == 0:
        return None
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
            
            if not content.strip():
                return None
            
            # Parse all timing values
            patterns = {
                'file_io': r'(?:off|ply)\s+read\s+file\s+I/O:\s+([\d.]+)\s+s',
                'read_vertices': r'(?:off|ply)\s+read\s+vertices:\s+([\d.]+)\s+s',
                'read_cells': r'(?:off|ply)\s+read\s+cells:\s+([\d.]+)\s+s',
                'build_ia_seq': r'build\s+IA\*\s+\(sequential\):\s+([\d.]+)\s+s',
                'build_ia_par': r'build\s+IA\*\s+\(parallel\):\s+([\d.]+)\s+s',
                'total_reading': r'(?:off|ply)\s+total\s+reading\s+time:\s+([\d.]+)\s+s',
                'gradient_encoding': r'Gradient\s+encoding\s+time:\s+([\d.]+)\s+s',
                'filtration': r'Filtration\s+computation\s+time:\s+([\d.]+)\s+s',
                'total_forman': r'Total\s+Forman\s+gradient\s+time:\s+([\d.]+)\s+s',
                'total_time': r'Total\s+time:\s+([\d.]+)\s+s',
                'total_reading_time': r'Total\s+reading\s+time:\s+([\d.]+)\s+s',
                'vertices': r'Vertices:\s+(\d+)',
                'top_simplexes': r'Total\s+top\s+simplexes:\s+(\d+)',
                'complex_dim': r'Complex\s+dimension:\s+(\d+)',
                'complex_vertices': r'Complex\s+vertices\s+#:\s+(\d+)',
                'complex_top_simplexes': r'Complex\s+top\s+simplices\s+#:\s+(\d+)',
                'reading_time': r'Reading\s+time\s+([\d.]+)\s+s',
                'filtration_time': r'Filtration\s+time\s+([\d.]+)\s+s',
                'forman_gradient_computed': r'Forman\s+gradient\s+computed\s+([\d.]+)\s+s'
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    try:
                        if key in ['vertices', 'top_simplexes', 'complex_vertices', 'complex_top_simplexes']:
                            data[key] = int(match.group(1))
                        else:
                            data[key] = float(match.group(1))
                    except ValueError:
                        pass
            
            # Also try to get build IA* time (without seq/par label)
            if 'build_ia_seq' not in data and 'build_ia_par' not in data:
                match = re.search(r'build\s+IA\*:\s+([\d.]+)\s+s', content, re.IGNORECASE)
                if match:
                    data['build_ia'] = float(match.group(1))
            
    except Exception as e:
        print(f"Error parsing {log_file}: {e}", file=sys.stderr)
        return None
    
    return data if data else None

def analyze_results(test_results_dir):
    """Analyze all test results in the directory."""
    
    print("="*100)
    print("TEST RESULTS ANALYSIS")
    print("="*100)
    print(f"Directory: {test_results_dir}\n")
    
    if not os.path.exists(test_results_dir):
        print(f"Error: Directory not found: {test_results_dir}")
        return
    
    # Collect all results
    results = defaultdict(dict)
    
    # Process .off files
    off_tests_dir = os.path.join(test_results_dir, "off_tests")
    if os.path.exists(off_tests_dir):
        print("Processing .off file tests...")
        for test_file_dir in os.listdir(off_tests_dir):
            test_file_path = os.path.join(off_tests_dir, test_file_dir)
            if os.path.isdir(test_file_path):
                print(f"  {test_file_dir}")
                for log_file in os.listdir(test_file_path):
                    if log_file.endswith('.log'):
                        # Handle both xx_tts and test_forman_gradient variants
                        variant = log_file.replace('_output.log', '')
                        log_path = os.path.join(test_file_path, log_file)
                        data = parse_timing_log(log_path)
                        if data:
                            results[test_file_dir][variant] = data
                        else:
                            file_size = os.path.getsize(log_path)
                            if file_size == 0:
                                print(f"    ⚠ {variant}: Empty log file")
                            else:
                                print(f"    ⚠ {variant}: Could not parse timing data")
    
    # Process .ply files
    ply_tests_dir = os.path.join(test_results_dir, "ply_tests")
    if os.path.exists(ply_tests_dir):
        print("\nProcessing .ply file tests...")
        for log_file in os.listdir(ply_tests_dir):
            if log_file.endswith('.log'):
                variant = log_file.replace('_output.log', '')
                log_path = os.path.join(ply_tests_dir, log_file)
                data = parse_timing_log(log_path)
                if data:
                    results['ply_file'][variant] = data
                else:
                    file_size = os.path.getsize(log_path)
                    if file_size == 0:
                        print(f"  ⚠ {variant}: Empty log file")
                    else:
                        print(f"  ⚠ {variant}: Could not parse timing data")
    
    if not results:
        print("\n⚠ No timing data found in log files!")
        print("\nPossible reasons:")
        print("1. Log files are empty (executables may need different arguments)")
        print("2. Timing format doesn't match expected patterns")
        print("3. Tests haven't been run yet")
        return
    
    # Generate summary tables
    print("\n" + "="*100)
    print("PERFORMANCE SUMMARY")
    print("="*100)
    
    for test_file, variants in sorted(results.items()):
        print(f"\n{test_file}")
        print("-" * 100)
        
        # Determine available metrics
        all_metrics = set()
        for variant_data in variants.values():
            all_metrics.update(variant_data.keys())
        
        # Filter to timing metrics
        timing_metrics = [m for m in ['build_ia_seq', 'build_ia_par', 'build_ia', 'total_reading', 
                                      'total_forman', 'total_time', 'reading_time', 
                                      'forman_gradient_computed'] if m in all_metrics]
        
        if timing_metrics:
            header = f"{'Variant':<20}"
            for metric in timing_metrics:
                header += f" {metric.replace('_', ' ').title():<15}"
            print(header)
            print("-" * 100)
            
            for variant, data in sorted(variants.items()):
                row = f"{variant:<20}"
                for metric in timing_metrics:
                    value = data.get(metric)
                    if value is not None:
                        row += f" {value:>14.3f}s"
                    else:
                        row += f" {'N/A':>15}"
                print(row)
        
        # File statistics
        if any('vertices' in data or 'complex_vertices' in data for data in variants.values()):
            print("\nFile Statistics:")
            print(f"{'Variant':<20} {'Vertices':<12} {'Top Simplexes':<15} {'Complex Dim':<12}")
            print("-" * 100)
            for variant, data in sorted(variants.items()):
                vertices = data.get('vertices') or data.get('complex_vertices')
                top_simplexes = data.get('top_simplexes') or data.get('complex_top_simplexes')
                dim = data.get('complex_dim')
                
                vertices_str = str(int(vertices)) if vertices else "N/A"
                top_str = str(int(top_simplexes)) if top_simplexes else "N/A"
                dim_str = str(int(dim)) if dim else "N/A"
                
                print(f"{variant:<20} {vertices_str:<12} {top_str:<15} {dim_str:<12}")
    
    # Calculate speedups if we have old and new data
    print("\n" + "="*100)
    print("SPEEDUP ANALYSIS (New vs Old)")
    print("="*100)
    
    for test_file, variants in sorted(results.items()):
        # Check for both xx_tts and test_forman_gradient variants
        for prefix in ['xx_tts', 'test_forman_gradient']:
            seq_old_key = f'{prefix}_seq_old'
            seq_new_key = f'{prefix}_seq_new'
            pa_old_key = f'{prefix}_pa_old'
            pa_new_key = f'{prefix}_pa_new'
            
            if seq_old_key in variants and seq_new_key in variants:
                print(f"\n{test_file} - {prefix} Sequential:")
                old_data = variants[seq_old_key]
                new_data = variants[seq_new_key]
                
                metrics_to_compare = ['build_ia_seq', 'build_ia_par', 'build_ia', 'total_reading', 
                                     'total_forman', 'total_time', 'reading_time']
                
                for metric in metrics_to_compare:
                    old_val = old_data.get(metric)
                    new_val = new_data.get(metric)
                    if old_val and new_val:
                        speedup = old_val / new_val
                        print(f"  {metric}: {old_val:.3f}s → {new_val:.3f}s (speedup: {speedup:.2f}x)")
            
            if pa_old_key in variants and pa_new_key in variants:
                print(f"\n{test_file} - {prefix} Parallel:")
                old_data = variants[pa_old_key]
                new_data = variants[pa_new_key]
                
                metrics_to_compare = ['build_ia_seq', 'build_ia_par', 'build_ia', 'total_reading', 
                                     'total_forman', 'total_time', 'reading_time']
                
                for metric in metrics_to_compare:
                    old_val = old_data.get(metric)
                    new_val = new_data.get(metric)
                    if old_val and new_val:
                        speedup = old_val / new_val
                        print(f"  {metric}: {old_val:.3f}s → {new_val:.3f}s (speedup: {speedup:.2f}x)")

def main():
    if len(sys.argv) < 2:
        print("Usage: python interpret_results.py <test_results_dir>")
        print("Example: python interpret_results.py cpp/cmake-build-debug/test_results_20251111_114634")
        sys.exit(1)
    
    test_results_dir = sys.argv[1]
    analyze_results(test_results_dir)

if __name__ == "__main__":
    main()

