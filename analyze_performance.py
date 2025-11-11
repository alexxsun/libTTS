#!/usr/bin/env python3
"""
Parse performance test results and generate analysis reports.

Usage:
    python analyze_performance.py [test_results_dir]
    
    Default: cpp/cmake-build-debug/test_results_*
"""

import sys
import os
import re
import glob
from pathlib import Path
from collections import defaultdict

def parse_timing_log(log_file):
    """Parse timing information from a log file."""
    data = {
        'file_io': None,
        'read_vertices': None,
        'read_cells': None,
        'build_ia': None,
        'total_reading': None,
        'gradient_encoding': None,
        'filtration': None,
        'total_forman': None,
        'total_time': None,
        'vertices': None,
        'top_simplexes': None,
        'complex_dim': None
    }
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
            
            # Parse timing values
            patterns = {
                'file_io': r'(?:off|ply)\s+read\s+file\s+I/O:\s+([\d.]+)\s+s',
                'read_vertices': r'(?:off|ply)\s+read\s+vertices:\s+([\d.]+)\s+s',
                'read_cells': r'(?:off|ply)\s+read\s+cells:\s+([\d.]+)\s+s',
                'build_ia': r'(?:off|ply)\s+build\s+IA\*\s+\((?:sequential|parallel)\):\s+([\d.]+)\s+s',
                'total_reading': r'(?:off|ply)\s+total\s+reading\s+time:\s+([\d.]+)\s+s',
                'gradient_encoding': r'Gradient\s+encoding\s+time:\s+([\d.]+)\s+s',
                'filtration': r'Filtration\s+computation\s+time:\s+([\d.]+)\s+s',
                'total_forman': r'Total\s+Forman\s+gradient\s+time:\s+([\d.]+)\s+s',
                'total_time': r'Total\s+time:\s+([\d.]+)\s+s',
                'vertices': r'Vertices:\s+(\d+)',
                'top_simplexes': r'Total\s+top\s+simplexes:\s+(\d+)',
                'complex_dim': r'Complex\s+dimension:\s+(\d+)'
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    try:
                        data[key] = float(match.group(1))
                    except ValueError:
                        pass
            
    except Exception as e:
        print(f"Error parsing {log_file}: {e}")
    
    return data

def find_test_results(base_dir):
    """Find all test result directories."""
    pattern = os.path.join(base_dir, "test_results_*")
    dirs = glob.glob(pattern)
    if not dirs:
        return None
    # Return most recent
    return max(dirs, key=os.path.getmtime)

def collect_performance_data(test_results_dir):
    """Collect all performance data from test results."""
    results = defaultdict(dict)
    
    # Find all log files
    off_tests_dir = os.path.join(test_results_dir, "off_tests")
    if os.path.exists(off_tests_dir):
        for test_file_dir in os.listdir(off_tests_dir):
            test_file_path = os.path.join(off_tests_dir, test_file_dir)
            if os.path.isdir(test_file_path):
                log_file = os.path.join(test_file_path, "test_forman_gradient.log")
                if os.path.exists(log_file):
                    data = parse_timing_log(log_file)
                    results[test_file_dir] = data
    
    return results

def generate_summary_table(results):
    """Generate a summary table of performance data."""
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY TABLE")
    print("="*80)
    
    if not results:
        print("No results found!")
        return
    
    # Table header
    header = f"{'File':<30} {'Vertices':<12} {'Build IA*':<12} {'Forman':<12} {'Total':<12}"
    print(header)
    print("-" * 80)
    
    for test_file, data in sorted(results.items()):
        vertices = int(data['vertices']) if data['vertices'] else 0
        build_ia = f"{data['build_ia']:.3f}" if data['build_ia'] else "N/A"
        forman = f"{data['total_forman']:.3f}" if data['total_forman'] else "N/A"
        total = f"{data['total_time']:.3f}" if data['total_time'] else "N/A"
        
        print(f"{test_file:<30} {vertices:<12} {build_ia:<12} {forman:<12} {total:<12}")

def generate_detailed_table(results):
    """Generate detailed timing breakdown."""
    print("\n" + "="*100)
    print("DETAILED TIMING BREAKDOWN")
    print("="*100)
    
    header = f"{'File':<30} {'File I/O':<10} {'Vertices':<10} {'Cells':<10} {'Build IA*':<12} {'Gradient':<10} {'Filtration':<12} {'Total':<10}"
    print(header)
    print("-" * 100)
    
    for test_file, data in sorted(results.items()):
        file_io = f"{data['file_io']:.3f}" if data['file_io'] else "N/A"
        vertices = f"{data['read_vertices']:.3f}" if data['read_vertices'] else "N/A"
        cells = f"{data['read_cells']:.3f}" if data['read_cells'] else "N/A"
        build_ia = f"{data['build_ia']:.3f}" if data['build_ia'] else "N/A"
        gradient = f"{data['gradient_encoding']:.3f}" if data['gradient_encoding'] else "N/A"
        filtration = f"{data['filtration']:.3f}" if data['filtration'] else "N/A"
        total = f"{data['total_time']:.3f}" if data['total_time'] else "N/A"
        
        print(f"{test_file:<30} {file_io:<10} {vertices:<10} {cells:<10} {build_ia:<12} {gradient:<10} {filtration:<12} {total:<10}")

def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "cpp/cmake-build-debug"
    
    print("="*80)
    print("Performance Analysis Tool")
    print("="*80)
    print(f"Searching in: {base_dir}")
    
    test_results_dir = find_test_results(base_dir)
    if not test_results_dir:
        print(f"Error: No test results found in {base_dir}")
        return 1
    
    print(f"Found test results: {os.path.basename(test_results_dir)}")
    
    results = collect_performance_data(test_results_dir)
    
    if not results:
        print("No performance data found in test results!")
        return 1
    
    print(f"\nFound data for {len(results)} test files")
    
    generate_summary_table(results)
    generate_detailed_table(results)
    
    # Calculate per-vertex metrics
    print("\n" + "="*80)
    print("PER-VERTEX METRICS")
    print("="*80)
    header = f"{'File':<30} {'Vertices':<12} {'Build/Vertex (ms)':<18} {'Forman/Vertex (ms)':<18}"
    print(header)
    print("-" * 80)
    
    for test_file, data in sorted(results.items()):
        vertices = data['vertices']
        if vertices and vertices > 0:
            build_ia = data['build_ia']
            forman = data['total_forman']
            
            build_per_vertex = (build_ia * 1000 / vertices) if build_ia else None
            forman_per_vertex = (forman * 1000 / vertices) if forman else None
            
            build_str = f"{build_per_vertex:.4f}" if build_per_vertex else "N/A"
            forman_str = f"{forman_per_vertex:.4f}" if forman_per_vertex else "N/A"
            
            print(f"{test_file:<30} {int(vertices):<12} {build_str:<18} {forman_str:<18}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

