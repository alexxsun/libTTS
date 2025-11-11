#!/usr/bin/env python3
"""
Create comprehensive performance report from test results.

Usage:
    python create_performance_report.py [base_dir]
    
    Default: cpp/cmake-build-debug
"""

import sys
import os
import re
import glob
from pathlib import Path
from collections import defaultdict

def parse_timing_from_log(log_file):
    """Parse timing information from a log file."""
    data = {}
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
            
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
                'complex_top_simplexes': r'Complex\s+top\s+simplices\s+#:\s+(\d+)'
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    try:
                        data[key] = float(match.group(1))
                    except ValueError:
                        pass
                # Also try integer match for counts
                if key in ['vertices', 'top_simplexes', 'complex_vertices', 'complex_top_simplexes']:
                    match = re.search(pattern, content, re.IGNORECASE)
                    if match:
                        try:
                            data[key] = int(match.group(1))
                        except ValueError:
                            pass
            
    except Exception as e:
        print(f"Error parsing {log_file}: {e}", file=sys.stderr)
    
    return data

def find_all_logs(base_dir):
    """Find all timing log files."""
    logs = {}
    
    # Check for direct log files
    for pattern in ['timing_*.log', '*_timing*.log']:
        for log_file in glob.glob(os.path.join(base_dir, pattern)):
            name = os.path.basename(log_file)
            logs[name] = log_file
    
    # Check test results directories
    for test_dir in glob.glob(os.path.join(base_dir, "test_results_*")):
        for root, dirs, files in os.walk(test_dir):
            for file in files:
                if file.endswith('.log'):
                    full_path = os.path.join(root, file)
                    rel_path = os.path.relpath(full_path, base_dir)
                    logs[rel_path] = full_path
    
    return logs

def generate_report(base_dir):
    """Generate comprehensive performance report."""
    print("="*100)
    print("PERFORMANCE TESTING REPORT")
    print("="*100)
    print(f"Base directory: {base_dir}\n")
    
    # Find all log files
    logs = find_all_logs(base_dir)
    
    if not logs:
        print("No log files found!")
        return
    
    print(f"Found {len(logs)} log file(s)\n")
    
    # Parse all logs
    all_data = {}
    for name, log_file in sorted(logs.items()):
        print(f"Parsing: {name}")
        data = parse_timing_from_log(log_file)
        if data:
            all_data[name] = data
    
    if not all_data:
        print("\nNo timing data found in log files!")
        return
    
    # Generate summary tables
    print("\n" + "="*100)
    print("TIMING SUMMARY")
    print("="*100)
    
    # Determine which metrics we have
    metrics = ['build_ia_seq', 'build_ia_par', 'total_reading', 'total_forman', 'total_time']
    available_metrics = []
    for metric in metrics:
        if any(metric in data for data in all_data.values()):
            available_metrics.append(metric)
    
    if available_metrics:
        header = f"{'Log File':<40}"
        for metric in available_metrics:
            header += f" {metric.replace('_', ' ').title():<15}"
        print(header)
        print("-" * 100)
        
        for name, data in sorted(all_data.items()):
            row = f"{name:<40}"
            for metric in available_metrics:
                value = data.get(metric)
                if value is not None:
                    row += f" {value:>14.3f}s"
                else:
                    row += f" {'N/A':>15}"
            print(row)
    
    # File statistics
    print("\n" + "="*100)
    print("FILE STATISTICS")
    print("="*100)
    
    stats_header = f"{'Log File':<40} {'Vertices':<12} {'Top Simplexes':<15} {'Complex Dim':<12}"
    print(stats_header)
    print("-" * 100)
    
    for name, data in sorted(all_data.items()):
        vertices = data.get('vertices') or data.get('complex_vertices')
        top_simplexes = data.get('top_simplexes') or data.get('complex_top_simplexes')
        dim = data.get('complex_dim')
        
        vertices_str = str(int(vertices)) if vertices else "N/A"
        top_str = str(int(top_simplexes)) if top_simplexes else "N/A"
        dim_str = str(int(dim)) if dim else "N/A"
        
        print(f"{name:<40} {vertices_str:<12} {top_str:<15} {dim_str:<12}")
    
    # Detailed breakdown
    print("\n" + "="*100)
    print("DETAILED TIMING BREAKDOWN")
    print("="*100)
    
    detail_header = f"{'Log File':<30} {'File I/O':<10} {'Vertices':<10} {'Cells':<10} {'Build IA*':<12} {'Gradient':<10} {'Filtration':<12} {'Total':<10}"
    print(detail_header)
    print("-" * 100)
    
    for name, data in sorted(all_data.items()):
        file_io = data.get('file_io')
        read_vertices = data.get('read_vertices')
        read_cells = data.get('read_cells')
        build_ia = data.get('build_ia_seq') or data.get('build_ia_par')
        gradient = data.get('gradient_encoding')
        filtration = data.get('filtration')
        total = data.get('total_time')
        
        name_short = name[:28] if len(name) > 28 else name
        
        file_io_str = f"{file_io:.3f}" if file_io else "N/A"
        vertices_str = f"{read_vertices:.3f}" if read_vertices else "N/A"
        cells_str = f"{read_cells:.3f}" if read_cells else "N/A"
        build_str = f"{build_ia:.3f}" if build_ia else "N/A"
        gradient_str = f"{gradient:.3f}" if gradient else "N/A"
        filtration_str = f"{filtration:.3f}" if filtration else "N/A"
        total_str = f"{total:.3f}" if total else "N/A"
        
        print(f"{name_short:<30} {file_io_str:<10} {vertices_str:<10} {cells_str:<10} {build_str:<12} {gradient_str:<10} {filtration_str:<12} {total_str:<10}")
    
    print("\n" + "="*100)
    print("Report generated successfully!")
    print("="*100)

if __name__ == "__main__":
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "cpp/cmake-build-debug"
    generate_report(base_dir)

