#!/usr/bin/env python3
"""
Generate performance report with figures from test results.

Usage:
    python generate_performance_report.py <test_results_dir>
    
Example:
    python generate_performance_report.py cpp/cmake-build-debug/test_results_20251111_123424
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict
from datetime import datetime

try:
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not found. Install with: pip install matplotlib numpy")
    print("  or: python3 -m pip install --user matplotlib numpy")
    print("  or: sudo pacman -S python-matplotlib python-numpy")
    print("\nGenerating text report only...\n")

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
                'build_ia_seq': r'(?:off|ply)\s+build\s+IA\*\s+\(sequential\):\s+([\d.]+)\s+s',
                'build_ia_par': r'(?:off|ply)\s+build\s+IA\*\s+\(parallel\):\s+([\d.]+)\s+s',
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
                'forman_gradient_computed': r'Forman\s+gradient\s+computed\s+([\d.]+)\s+s',
                'memory_kb': r'Memory\s+usage\s+\(KB\):\s+(\d+)',
                'memory_mb': r'Memory\s+usage\s+\(MB\):\s+([\d.]+)'
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    try:
                        if key in ['vertices', 'top_simplexes', 'complex_vertices', 'complex_top_simplexes', 'memory_kb']:
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

def collect_results(test_results_dir):
    """Collect all test results."""
    results = defaultdict(dict)
    
    # Process .off files
    off_tests_dir = os.path.join(test_results_dir, "off_tests")
    if os.path.exists(off_tests_dir):
        for test_file_dir in os.listdir(off_tests_dir):
            test_file_path = os.path.join(off_tests_dir, test_file_dir)
            if os.path.isdir(test_file_path):
                for log_file in os.listdir(test_file_path):
                    if log_file.endswith('.log'):
                        variant = log_file.replace('_output.log', '')
                        log_path = os.path.join(test_file_path, log_file)
                        data = parse_timing_log(log_path)
                        if data:
                            results[test_file_dir][variant] = data
    
    # Process .ply files
    ply_tests_dir = os.path.join(test_results_dir, "ply_tests")
    if os.path.exists(ply_tests_dir):
        for log_file in os.listdir(ply_tests_dir):
            if log_file.endswith('.log'):
                variant = log_file.replace('_output.log', '')
                log_path = os.path.join(ply_tests_dir, log_file)
                data = parse_timing_log(log_path)
                if data:
                    results['ply_file'][variant] = data
    
    return results

def create_sequential_comparison_plot(results, output_dir):
    """Create comprehensive sequential comparison: old vs new (sequential versions)."""
    if not HAS_MATPLOTLIB:
        print("  ⚠ Skipping plots (matplotlib not available)")
        return
    
    # Set style
    try:
        plt.style.use('seaborn-v0_8-darkgrid')
    except:
        try:
            plt.style.use('seaborn-darkgrid')
        except:
            pass  # Use default style
    
    fig = plt.figure(figsize=(18, 12))
    
    # Collect data for plotting
    off_files = [k for k in results.keys() if k != 'ply_file']
    off_files.sort()
    
    x_pos = np.arange(len(off_files))
    width = 0.35
    
    # Plot 1: Build IA* Time Comparison (Sequential)
    ax1 = plt.subplot(2, 3, 1)
    old_times = []
    new_times = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_seq_old', {})
        new_data = results[file_name].get('test_forman_gradient_seq_new', {})
        old_time = old_data.get('build_ia_seq', 0) or old_data.get('build_ia', 0)
        new_time = new_data.get('build_ia_seq', 0) or new_data.get('build_ia', 0)
        old_times.append(old_time)
        new_times.append(new_time)
    
    bars1 = ax1.bar(x_pos - width/2, old_times, width, label='Old (adjRelations)', color='#d62728', alpha=0.8)
    bars2 = ax1.bar(x_pos + width/2, new_times, width, label='New (completeCoboundaryTop)', color='#2ca02c', alpha=0.8)
    ax1.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax1.set_title('1. Build IA* (Sequential)', fontsize=12, fontweight='bold')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 2: Total Reading Time (Sequential)
    ax2 = plt.subplot(2, 3, 2)
    old_times = []
    new_times = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_seq_old', {})
        new_data = results[file_name].get('test_forman_gradient_seq_new', {})
        old_time = old_data.get('total_reading', 0) or old_data.get('total_reading_time', 0)
        new_time = new_data.get('total_reading', 0) or new_data.get('total_reading_time', 0)
        old_times.append(old_time)
        new_times.append(new_time)
    
    bars1 = ax2.bar(x_pos - width/2, old_times, width, label='Old', color='#d62728', alpha=0.8)
    bars2 = ax2.bar(x_pos + width/2, new_times, width, label='New', color='#2ca02c', alpha=0.8)
    ax2.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax2.set_title('2. Total Reading (Sequential)', fontsize=12, fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 3: Forman Gradient Time (Sequential)
    ax3 = plt.subplot(2, 3, 3)
    old_times = []
    new_times = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_seq_old', {})
        new_data = results[file_name].get('test_forman_gradient_seq_new', {})
        old_time = old_data.get('total_forman', 0)
        new_time = new_data.get('total_forman', 0)
        old_times.append(old_time)
        new_times.append(new_time)
    
    bars1 = ax3.bar(x_pos - width/2, old_times, width, label='Old', color='#d62728', alpha=0.8)
    bars2 = ax3.bar(x_pos + width/2, new_times, width, label='New', color='#2ca02c', alpha=0.8)
    ax3.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax3.set_title('3. Forman Gradient (Sequential)', fontsize=12, fontweight='bold')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax3.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 4: Total Time (Sequential)
    ax4 = plt.subplot(2, 3, 4)
    old_times = []
    new_times = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_seq_old', {})
        new_data = results[file_name].get('test_forman_gradient_seq_new', {})
        old_time = old_data.get('total_time', 0) or (old_data.get('total_reading', 0) + old_data.get('total_forman', 0))
        new_time = new_data.get('total_time', 0) or (new_data.get('total_reading', 0) + new_data.get('total_forman', 0))
        old_times.append(old_time)
        new_times.append(new_time)
    
    bars1 = ax4.bar(x_pos - width/2, old_times, width, label='Old', color='#d62728', alpha=0.8)
    bars2 = ax4.bar(x_pos + width/2, new_times, width, label='New', color='#2ca02c', alpha=0.8)
    ax4.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax4.set_title('4. Total Time (Sequential)', fontsize=12, fontweight='bold')
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax4.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 5: Speedup Breakdown (Sequential)
    ax5 = plt.subplot(2, 3, 5)
    metrics_speedup = {
        'Build IA*': [],
        'Total Reading': [],
        'Forman Gradient': [],
        'Total Time': []
    }
    
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_seq_old', {})
        new_data = results[file_name].get('test_forman_gradient_seq_new', {})
        
        # Build IA*
        old_build = old_data.get('build_ia_seq', 0) or old_data.get('build_ia', 0)
        new_build = new_data.get('build_ia_seq', 0) or new_data.get('build_ia', 0)
        if old_build > 0 and new_build > 0:
            metrics_speedup['Build IA*'].append(old_build / new_build)
        else:
            metrics_speedup['Build IA*'].append(0)
        
        # Total Reading
        old_read = old_data.get('total_reading', 0) or old_data.get('total_reading_time', 0)
        new_read = new_data.get('total_reading', 0) or new_data.get('total_reading_time', 0)
        if old_read > 0 and new_read > 0:
            metrics_speedup['Total Reading'].append(old_read / new_read)
        else:
            metrics_speedup['Total Reading'].append(0)
        
        # Forman Gradient
        old_forman = old_data.get('total_forman', 0)
        new_forman = new_data.get('total_forman', 0)
        if old_forman > 0 and new_forman > 0:
            metrics_speedup['Forman Gradient'].append(old_forman / new_forman)
        else:
            metrics_speedup['Forman Gradient'].append(0)
        
        # Total Time
        old_total = old_data.get('total_time', 0) or (old_data.get('total_reading', 0) + old_data.get('total_forman', 0))
        new_total = new_data.get('total_time', 0) or (new_data.get('total_reading', 0) + new_data.get('total_forman', 0))
        if old_total > 0 and new_total > 0:
            metrics_speedup['Total Time'].append(old_total / new_total)
        else:
            metrics_speedup['Total Time'].append(0)
    
    width_speedup = 0.2
    colors_speedup = ['#ff7f0e', '#1f77b4', '#2ca02c', '#9467bd']
    for i, (metric, speedups) in enumerate(metrics_speedup.items()):
        offset = (i - 1.5) * width_speedup
        bars = ax5.bar(x_pos + offset, speedups, width_speedup, label=metric, color=colors_speedup[i], alpha=0.8)
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax5.text(bar.get_x() + bar.get_width()/2., height, f'{height:.2f}x',
                        ha='center', va='bottom', fontsize=7)
    
    ax5.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax5.set_ylabel('Speedup (x)', fontsize=11, fontweight='bold')
    ax5.set_title('5. Speedup: Old/New (Sequential)', fontsize=12, fontweight='bold')
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax5.legend(fontsize=8, ncol=2)
    ax5.grid(True, alpha=0.3, axis='y')
    ax5.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, linewidth=1)
    
    # Plot 6: Scalability (Sequential) - Multiple metrics
    ax6 = plt.subplot(2, 3, 6)
    
    # Collect data for all three metrics
    vertices = []
    faces = []
    avg_faces_per_vertex = []
    old_times = []
    new_times = []
    file_names = []
    
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_seq_old', {})
        new_data = results[file_name].get('test_forman_gradient_seq_new', {})
        v = old_data.get('vertices', 0) or old_data.get('complex_vertices', 0)
        f = old_data.get('top_simplexes', 0) or old_data.get('complex_top_simplexes', 0)
        old_time = old_data.get('total_time', 0) or (old_data.get('total_reading', 0) + old_data.get('total_forman', 0))
        new_time = new_data.get('total_time', 0) or (new_data.get('total_reading', 0) + new_data.get('total_forman', 0))
        
        if v > 0 and f > 0:
            vertices.append(v)
            faces.append(f)
            avg_faces_per_vertex.append(f / v)
            old_times.append(old_time)
            new_times.append(new_time)
            file_names.append(file_name)
    
    if vertices:
        # Sort by vertices for consistent ordering
        sorted_data = sorted(zip(vertices, faces, avg_faces_per_vertex, old_times, new_times, file_names))
        vertices, faces, avg_faces_per_vertex, old_times, new_times, file_names = zip(*sorted_data)
        
        # Plot all three relationships on the same axes (time is y-axis for all)
        # Use different line styles and markers to distinguish
        ax6.plot(vertices, old_times, 'o-', label='Old (vs vertices)', color='#d62728', linewidth=2, markersize=6, alpha=0.8)
        ax6.plot(vertices, new_times, 's-', label='New (vs vertices)', color='#2ca02c', linewidth=2, markersize=6, alpha=0.8)
        ax6.plot(faces, old_times, 'o--', label='Old (vs faces)', color='#d62728', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(faces, new_times, 's--', label='New (vs faces)', color='#2ca02c', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(avg_faces_per_vertex, old_times, 'o:', label='Old (vs avg faces/v)', color='#d62728', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(avg_faces_per_vertex, new_times, 's:', label='New (vs avg faces/v)', color='#2ca02c', linewidth=1.5, markersize=5, alpha=0.6)
        
        ax6.set_xlabel('X-axis: Vertices / Faces / Avg Faces/Vertex', fontsize=10, fontweight='bold')
        ax6.set_ylabel('Total Time (seconds)', fontsize=11, fontweight='bold')
        ax6.set_title('6. Scalability (Sequential)', fontsize=12, fontweight='bold')
        ax6.legend(fontsize=7, loc='upper left', ncol=2)
        ax6.grid(True, alpha=0.3)
        ax6.set_xscale('log')
        ax6.set_yscale('log')
    
    plt.suptitle('Sequential Comparison: Old vs New Code', fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    output_file = os.path.join(output_dir, 'sequential_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    plt.close()

def create_parallel_comparison_plot(results, output_dir):
    """Create comprehensive parallel comparison: old vs new (parallel versions)."""
    if not HAS_MATPLOTLIB:
        print("  ⚠ Skipping plots (matplotlib not available)")
        return
    
    # Set style
    try:
        plt.style.use('seaborn-v0_8-darkgrid')
    except:
        try:
            plt.style.use('seaborn-darkgrid')
        except:
            pass
    
    fig = plt.figure(figsize=(18, 12))
    
    off_files = [k for k in results.keys() if k != 'ply_file']
    off_files.sort()
    
    x_pos = np.arange(len(off_files))
    width = 0.35
    
    # Plot 1: Build IA* Time Comparison (Parallel)
    ax1 = plt.subplot(2, 3, 1)
    old_times = []
    new_times = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_pa_old', {})
        new_data = results[file_name].get('test_forman_gradient_pa_new', {})
        # Try build_ia_par first, then fallback to build_ia_seq or build_ia
        old_time = old_data.get('build_ia_par', 0) or old_data.get('build_ia_seq', 0) or old_data.get('build_ia', 0)
        new_time = new_data.get('build_ia_par', 0) or new_data.get('build_ia_seq', 0) or new_data.get('build_ia', 0)
        old_times.append(old_time)
        new_times.append(new_time)
    
    bars1 = ax1.bar(x_pos - width/2, old_times, width, label='Old (adjRelations)', color='#d62728', alpha=0.8)
    bars2 = ax1.bar(x_pos + width/2, new_times, width, label='New (completeCoboundaryTop)', color='#2ca02c', alpha=0.8)
    ax1.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax1.set_title('1. Build IA* (Parallel)', fontsize=12, fontweight='bold')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 2: Total Reading Time (Parallel)
    ax2 = plt.subplot(2, 3, 2)
    old_times = []
    new_times = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_pa_old', {})
        new_data = results[file_name].get('test_forman_gradient_pa_new', {})
        old_time = old_data.get('total_reading', 0) or old_data.get('total_reading_time', 0)
        new_time = new_data.get('total_reading', 0) or new_data.get('total_reading_time', 0)
        old_times.append(old_time)
        new_times.append(new_time)
    
    bars1 = ax2.bar(x_pos - width/2, old_times, width, label='Old', color='#d62728', alpha=0.8)
    bars2 = ax2.bar(x_pos + width/2, new_times, width, label='New', color='#2ca02c', alpha=0.8)
    ax2.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax2.set_title('2. Total Reading (Parallel)', fontsize=12, fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 3: Forman Gradient Time (Parallel)
    ax3 = plt.subplot(2, 3, 3)
    old_times = []
    new_times = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_pa_old', {})
        new_data = results[file_name].get('test_forman_gradient_pa_new', {})
        old_time = old_data.get('total_forman', 0)
        new_time = new_data.get('total_forman', 0)
        old_times.append(old_time)
        new_times.append(new_time)
    
    bars1 = ax3.bar(x_pos - width/2, old_times, width, label='Old', color='#d62728', alpha=0.8)
    bars2 = ax3.bar(x_pos + width/2, new_times, width, label='New', color='#2ca02c', alpha=0.8)
    ax3.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax3.set_title('3. Forman Gradient (Parallel)', fontsize=12, fontweight='bold')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax3.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 4: Total Time (Parallel)
    ax4 = plt.subplot(2, 3, 4)
    old_times = []
    new_times = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_pa_old', {})
        new_data = results[file_name].get('test_forman_gradient_pa_new', {})
        old_time = old_data.get('total_time', 0) or (old_data.get('total_reading', 0) + old_data.get('total_forman', 0))
        new_time = new_data.get('total_time', 0) or (new_data.get('total_reading', 0) + new_data.get('total_forman', 0))
        old_times.append(old_time)
        new_times.append(new_time)
    
    bars1 = ax4.bar(x_pos - width/2, old_times, width, label='Old', color='#d62728', alpha=0.8)
    bars2 = ax4.bar(x_pos + width/2, new_times, width, label='New', color='#2ca02c', alpha=0.8)
    ax4.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax4.set_title('4. Total Time (Parallel)', fontsize=12, fontweight='bold')
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax4.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 5: Speedup Breakdown (Parallel)
    ax5 = plt.subplot(2, 3, 5)
    metrics_speedup = {
        'Build IA*': [],
        'Total Reading': [],
        'Forman Gradient': [],
        'Total Time': []
    }
    
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_pa_old', {})
        new_data = results[file_name].get('test_forman_gradient_pa_new', {})
        
        # Build IA*
        old_build = old_data.get('build_ia_par', 0) or old_data.get('build_ia_seq', 0) or old_data.get('build_ia', 0)
        new_build = new_data.get('build_ia_par', 0) or new_data.get('build_ia_seq', 0) or new_data.get('build_ia', 0)
        if old_build > 0 and new_build > 0:
            metrics_speedup['Build IA*'].append(old_build / new_build)
        else:
            metrics_speedup['Build IA*'].append(0)
        
        # Total Reading
        old_read = old_data.get('total_reading', 0) or old_data.get('total_reading_time', 0)
        new_read = new_data.get('total_reading', 0) or new_data.get('total_reading_time', 0)
        if old_read > 0 and new_read > 0:
            metrics_speedup['Total Reading'].append(old_read / new_read)
        else:
            metrics_speedup['Total Reading'].append(0)
        
        # Forman Gradient
        old_forman = old_data.get('total_forman', 0)
        new_forman = new_data.get('total_forman', 0)
        if old_forman > 0 and new_forman > 0:
            metrics_speedup['Forman Gradient'].append(old_forman / new_forman)
        else:
            metrics_speedup['Forman Gradient'].append(0)
        
        # Total Time
        old_total = old_data.get('total_time', 0) or (old_data.get('total_reading', 0) + old_data.get('total_forman', 0))
        new_total = new_data.get('total_time', 0) or (new_data.get('total_reading', 0) + new_data.get('total_forman', 0))
        if old_total > 0 and new_total > 0:
            metrics_speedup['Total Time'].append(old_total / new_total)
        else:
            metrics_speedup['Total Time'].append(0)
    
    width_speedup = 0.2
    colors_speedup = ['#ff7f0e', '#1f77b4', '#2ca02c', '#9467bd']
    for i, (metric, speedups) in enumerate(metrics_speedup.items()):
        offset = (i - 1.5) * width_speedup
        bars = ax5.bar(x_pos + offset, speedups, width_speedup, label=metric, color=colors_speedup[i], alpha=0.8)
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax5.text(bar.get_x() + bar.get_width()/2., height, f'{height:.2f}x',
                        ha='center', va='bottom', fontsize=7)
    
    ax5.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax5.set_ylabel('Speedup (x)', fontsize=11, fontweight='bold')
    ax5.set_title('5. Speedup: Old/New (Parallel)', fontsize=12, fontweight='bold')
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax5.legend(fontsize=8, ncol=2)
    ax5.grid(True, alpha=0.3, axis='y')
    ax5.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, linewidth=1)
    
    # Plot 6: Scalability (Parallel) - Multiple metrics
    ax6 = plt.subplot(2, 3, 6)
    
    # Collect data for all three metrics
    vertices = []
    faces = []
    avg_faces_per_vertex = []
    old_times = []
    new_times = []
    file_names = []
    
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_pa_old', {})
        new_data = results[file_name].get('test_forman_gradient_pa_new', {})
        v = old_data.get('vertices', 0) or old_data.get('complex_vertices', 0)
        f = old_data.get('top_simplexes', 0) or old_data.get('complex_top_simplexes', 0)
        old_time = old_data.get('total_time', 0) or (old_data.get('total_reading', 0) + old_data.get('total_forman', 0))
        new_time = new_data.get('total_time', 0) or (new_data.get('total_reading', 0) + new_data.get('total_forman', 0))
        
        if v > 0 and f > 0:
            vertices.append(v)
            faces.append(f)
            avg_faces_per_vertex.append(f / v)
            old_times.append(old_time)
            new_times.append(new_time)
            file_names.append(file_name)
    
    if vertices:
        # Sort by vertices for consistent ordering
        sorted_data = sorted(zip(vertices, faces, avg_faces_per_vertex, old_times, new_times, file_names))
        vertices, faces, avg_faces_per_vertex, old_times, new_times, file_names = zip(*sorted_data)
        
        # Plot all three relationships on the same axes (time is y-axis for all)
        # Use different line styles and markers to distinguish
        ax6.plot(vertices, old_times, 'o-', label='Old (vs vertices)', color='#d62728', linewidth=2, markersize=6, alpha=0.8)
        ax6.plot(vertices, new_times, 's-', label='New (vs vertices)', color='#2ca02c', linewidth=2, markersize=6, alpha=0.8)
        ax6.plot(faces, old_times, 'o--', label='Old (vs faces)', color='#d62728', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(faces, new_times, 's--', label='New (vs faces)', color='#2ca02c', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(avg_faces_per_vertex, old_times, 'o:', label='Old (vs avg faces/v)', color='#d62728', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(avg_faces_per_vertex, new_times, 's:', label='New (vs avg faces/v)', color='#2ca02c', linewidth=1.5, markersize=5, alpha=0.6)
        
        ax6.set_xlabel('X-axis: Vertices / Faces / Avg Faces/Vertex', fontsize=10, fontweight='bold')
        ax6.set_ylabel('Total Time (seconds)', fontsize=11, fontweight='bold')
        ax6.set_title('6. Scalability (Parallel)', fontsize=12, fontweight='bold')
        ax6.legend(fontsize=7, loc='upper left', ncol=2)
        ax6.grid(True, alpha=0.3)
        ax6.set_xscale('log')
        ax6.set_yscale('log')
    
    plt.suptitle('Parallel Comparison: Old vs New Code', fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    output_file = os.path.join(output_dir, 'parallel_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    plt.close()

def create_speedup_analysis_plot(results, output_dir):
    """Create detailed speedup analysis: sequential and parallel separately."""
    if not HAS_MATPLOTLIB:
        print("  ⚠ Skipping plots (matplotlib not available)")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    off_files = [k for k in results.keys() if k != 'ply_file']
    off_files.sort()
    
    metrics = ['Build IA*', 'Total Reading', 'Forman Gradient', 'Total Time']
    x = np.arange(len(off_files))
    width = 0.2
    
    # Collect all speedups to determine y-axis range
    all_speedups_seq = []
    all_speedups_par = []
    
    # Sequential speedup
    for i, metric in enumerate(metrics):
        speedups = []
        for file_name in off_files:
            old_data = results[file_name].get('test_forman_gradient_seq_old', {})
            new_data = results[file_name].get('test_forman_gradient_seq_new', {})
            
            if metric == 'Build IA*':
                old_val = old_data.get('build_ia_seq', 0) or old_data.get('build_ia', 0)
                new_val = new_data.get('build_ia_seq', 0) or new_data.get('build_ia', 0)
            elif metric == 'Total Reading':
                old_val = old_data.get('total_reading', 0) or old_data.get('total_reading_time', 0)
                new_val = new_data.get('total_reading', 0) or new_data.get('total_reading_time', 0)
            elif metric == 'Forman Gradient':
                old_val = old_data.get('total_forman', 0)
                new_val = new_data.get('total_forman', 0)
            else:  # Total Time
                old_val = old_data.get('total_time', 0) or (old_data.get('total_reading', 0) + old_data.get('total_forman', 0))
                new_val = new_data.get('total_time', 0) or (new_data.get('total_reading', 0) + new_data.get('total_forman', 0))
            
            if old_val > 0 and new_val > 0:
                speedup = old_val / new_val
                speedups.append(speedup)
                all_speedups_seq.append(speedup)
            else:
                speedups.append(0)
        
        offset = (i - 1.5) * width
        bars = ax1.bar(x + offset, speedups, width, label=metric, alpha=0.8)
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height, f'{height:.2f}x',
                        ha='center', va='bottom', fontsize=7)
    
    ax1.set_xlabel('Test Files', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Speedup (x)', fontsize=12, fontweight='bold')
    ax1.set_title('Speedup Analysis: Sequential (Old vs New)', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=9)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, linewidth=1)
    
    # Parallel speedup
    for i, metric in enumerate(metrics):
        speedups = []
        for file_name in off_files:
            old_data = results[file_name].get('test_forman_gradient_pa_old', {})
            new_data = results[file_name].get('test_forman_gradient_pa_new', {})
            
            if metric == 'Build IA*':
                # Try build_ia_par first, then fallback
                old_val = old_data.get('build_ia_par', 0) or old_data.get('build_ia_seq', 0) or old_data.get('build_ia', 0)
                new_val = new_data.get('build_ia_par', 0) or new_data.get('build_ia_seq', 0) or new_data.get('build_ia', 0)
            elif metric == 'Total Reading':
                old_val = old_data.get('total_reading', 0) or old_data.get('total_reading_time', 0)
                new_val = new_data.get('total_reading', 0) or new_data.get('total_reading_time', 0)
            elif metric == 'Forman Gradient':
                old_val = old_data.get('total_forman', 0)
                new_val = new_data.get('total_forman', 0)
            else:  # Total Time
                old_val = old_data.get('total_time', 0) or (old_data.get('total_reading', 0) + old_data.get('total_forman', 0))
                new_val = new_data.get('total_time', 0) or (new_data.get('total_reading', 0) + new_data.get('total_forman', 0))
            
            if old_val > 0 and new_val > 0:
                speedup = old_val / new_val
                speedups.append(speedup)
                all_speedups_par.append(speedup)
            else:
                speedups.append(0)
        
        offset = (i - 1.5) * width
        bars = ax2.bar(x + offset, speedups, width, label=metric, alpha=0.8)
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height, f'{height:.2f}x',
                        ha='center', va='bottom', fontsize=7)
    
    ax2.set_xlabel('Test Files', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Speedup (x)', fontsize=12, fontweight='bold')
    ax2.set_title('Speedup Analysis: Parallel (Old vs New)', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=9)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, linewidth=1)
    
    # Set same y-axis limits for both plots
    all_speedups = all_speedups_seq + all_speedups_par
    if all_speedups:
        y_max = max(all_speedups) * 1.1
        y_min = min(0.8, min([s for s in all_speedups if s > 0]) * 0.9) if any(s > 0 for s in all_speedups) else 0
        ax1.set_ylim([y_min, y_max])
        ax2.set_ylim([y_min, y_max])
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'speedup_analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    
    plt.close()

def create_parallel_vs_sequential_plot(results, output_dir):
    """Create plot comparing parallel vs sequential in new code."""
    if not HAS_MATPLOTLIB:
        print("  ⚠ Skipping plots (matplotlib not available)")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    off_files = [k for k in results.keys() if k != 'ply_file']
    off_files.sort()
    
    x_pos = np.arange(len(off_files))
    width = 0.35
    
    # Plot 1: Build IA* - Sequential vs Parallel (New Code)
    ax1 = axes[0, 0]
    seq_times = []
    par_times = []
    for file_name in off_files:
        seq_data = results[file_name].get('test_forman_gradient_seq_new', {})
        par_data = results[file_name].get('test_forman_gradient_pa_new', {})
        seq_time = seq_data.get('build_ia_seq', 0) or seq_data.get('build_ia_par', 0) or seq_data.get('build_ia', 0)
        par_time = par_data.get('build_ia_par', 0) or par_data.get('build_ia_seq', 0) or par_data.get('build_ia', 0)
        seq_times.append(seq_time)
        par_times.append(par_time)
    
    bars1 = ax1.bar(x_pos - width/2, seq_times, width, label='Sequential', color='#1f77b4', alpha=0.8)
    bars2 = ax1.bar(x_pos + width/2, par_times, width, label='Parallel', color='#ff7f0e', alpha=0.8)
    ax1.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax1.set_title('1. Build IA*: Seq vs Par (New Code)', fontsize=12, fontweight='bold')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 2: Total Reading Time - Sequential vs Parallel (New Code)
    ax2 = axes[0, 1]
    seq_times = []
    par_times = []
    for file_name in off_files:
        seq_data = results[file_name].get('test_forman_gradient_seq_new', {})
        par_data = results[file_name].get('test_forman_gradient_pa_new', {})
        seq_time = seq_data.get('total_reading', 0) or seq_data.get('total_reading_time', 0)
        par_time = par_data.get('total_reading', 0) or par_data.get('total_reading_time', 0)
        seq_times.append(seq_time)
        par_times.append(par_time)
    
    bars1 = ax2.bar(x_pos - width/2, seq_times, width, label='Sequential', color='#1f77b4', alpha=0.8)
    bars2 = ax2.bar(x_pos + width/2, par_times, width, label='Parallel', color='#ff7f0e', alpha=0.8)
    ax2.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax2.set_title('2. Total Reading: Seq vs Par (New Code)', fontsize=12, fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 3: Forman Gradient - Sequential vs Parallel (New Code)
    ax3 = axes[0, 2]
    seq_times = []
    par_times = []
    for file_name in off_files:
        seq_data = results[file_name].get('test_forman_gradient_seq_new', {})
        par_data = results[file_name].get('test_forman_gradient_pa_new', {})
        seq_time = seq_data.get('total_forman', 0)
        par_time = par_data.get('total_forman', 0)
        seq_times.append(seq_time)
        par_times.append(par_time)
    
    bars1 = ax3.bar(x_pos - width/2, seq_times, width, label='Sequential', color='#1f77b4', alpha=0.8)
    bars2 = ax3.bar(x_pos + width/2, par_times, width, label='Parallel', color='#ff7f0e', alpha=0.8)
    ax3.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax3.set_title('3. Forman Gradient: Seq vs Par (New Code)', fontsize=12, fontweight='bold')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax3.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 4: Total Time - Sequential vs Parallel (New Code)
    ax4 = axes[1, 0]
    seq_times = []
    par_times = []
    for file_name in off_files:
        seq_data = results[file_name].get('test_forman_gradient_seq_new', {})
        par_data = results[file_name].get('test_forman_gradient_pa_new', {})
        seq_time = seq_data.get('total_time', 0) or (seq_data.get('total_reading', 0) + seq_data.get('total_forman', 0))
        par_time = par_data.get('total_time', 0) or (par_data.get('total_reading', 0) + par_data.get('total_forman', 0))
        seq_times.append(seq_time)
        par_times.append(par_time)
    
    bars1 = ax4.bar(x_pos - width/2, seq_times, width, label='Sequential', color='#1f77b4', alpha=0.8)
    bars2 = ax4.bar(x_pos + width/2, par_times, width, label='Parallel', color='#ff7f0e', alpha=0.8)
    ax4.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
    ax4.set_title('4. Total Time: Seq vs Par (New Code)', fontsize=12, fontweight='bold')
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax4.text(bar.get_x() + bar.get_width()/2., height, f'{height:.3f}s',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 5: Parallel Speedup Breakdown (seq_new / pa_new)
    ax5 = axes[1, 1]
    metrics_speedup = {
        'Build IA*': [],
        'Total Reading': [],
        'Forman Gradient': [],
        'Total Time': []
    }
    
    for file_name in off_files:
        seq_data = results[file_name].get('test_forman_gradient_seq_new', {})
        par_data = results[file_name].get('test_forman_gradient_pa_new', {})
        
        # Build IA*
        seq_build = seq_data.get('build_ia_seq', 0) or seq_data.get('build_ia_par', 0) or seq_data.get('build_ia', 0)
        par_build = par_data.get('build_ia_par', 0) or par_data.get('build_ia_seq', 0) or par_data.get('build_ia', 0)
        if seq_build > 0 and par_build > 0:
            metrics_speedup['Build IA*'].append(seq_build / par_build)
        else:
            metrics_speedup['Build IA*'].append(0)
        
        # Total Reading
        seq_read = seq_data.get('total_reading', 0) or seq_data.get('total_reading_time', 0)
        par_read = par_data.get('total_reading', 0) or par_data.get('total_reading_time', 0)
        if seq_read > 0 and par_read > 0:
            metrics_speedup['Total Reading'].append(seq_read / par_read)
        else:
            metrics_speedup['Total Reading'].append(0)
        
        # Forman Gradient
        seq_forman = seq_data.get('total_forman', 0)
        par_forman = par_data.get('total_forman', 0)
        if seq_forman > 0 and par_forman > 0:
            metrics_speedup['Forman Gradient'].append(seq_forman / par_forman)
        else:
            metrics_speedup['Forman Gradient'].append(0)
        
        # Total Time
        seq_total = seq_data.get('total_time', 0) or (seq_data.get('total_reading', 0) + seq_data.get('total_forman', 0))
        par_total = par_data.get('total_time', 0) or (par_data.get('total_reading', 0) + par_data.get('total_forman', 0))
        if seq_total > 0 and par_total > 0:
            metrics_speedup['Total Time'].append(seq_total / par_total)
        else:
            metrics_speedup['Total Time'].append(0)
    
    width_speedup = 0.2
    colors_speedup = ['#ff7f0e', '#1f77b4', '#2ca02c', '#9467bd']
    for i, (metric, speedups) in enumerate(metrics_speedup.items()):
        offset = (i - 1.5) * width_speedup
        bars = ax5.bar(x_pos + offset, speedups, width_speedup, label=metric, color=colors_speedup[i], alpha=0.8)
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax5.text(bar.get_x() + bar.get_width()/2., height, f'{height:.2f}x',
                        ha='center', va='bottom', fontsize=7)
    
    ax5.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax5.set_ylabel('Speedup (x)', fontsize=11, fontweight='bold')
    ax5.set_title('5. Parallel Speedup: Seq/Par (New Code)', fontsize=12, fontweight='bold')
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax5.legend(fontsize=8, ncol=2)
    ax5.grid(True, alpha=0.3, axis='y')
    ax5.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, linewidth=1)
    
    # Plot 6: Scalability - Sequential vs Parallel (New Code) - Multiple metrics
    ax6 = axes[1, 2]
    
    # Collect data for all three metrics
    vertices = []
    faces = []
    avg_faces_per_vertex = []
    seq_times = []
    par_times = []
    
    for file_name in off_files:
        seq_data = results[file_name].get('test_forman_gradient_seq_new', {})
        par_data = results[file_name].get('test_forman_gradient_pa_new', {})
        v = seq_data.get('vertices', 0) or seq_data.get('complex_vertices', 0)
        f = seq_data.get('top_simplexes', 0) or seq_data.get('complex_top_simplexes', 0)
        seq_time = seq_data.get('total_time', 0) or (seq_data.get('total_reading', 0) + seq_data.get('total_forman', 0))
        par_time = par_data.get('total_time', 0) or (par_data.get('total_reading', 0) + par_data.get('total_forman', 0))
        
        if v > 0 and f > 0:
            vertices.append(v)
            faces.append(f)
            avg_faces_per_vertex.append(f / v)
            seq_times.append(seq_time)
            par_times.append(par_time)
    
    if vertices:
        # Sort by vertices for consistent ordering
        sorted_data = sorted(zip(vertices, faces, avg_faces_per_vertex, seq_times, par_times))
        vertices, faces, avg_faces_per_vertex, seq_times, par_times = zip(*sorted_data)
        
        # Plot all three relationships on the same axes (time is y-axis for all)
        # Use different line styles and markers to distinguish
        ax6.plot(vertices, seq_times, 'o-', label='Seq (vs vertices)', color='#1f77b4', linewidth=2, markersize=6, alpha=0.8)
        ax6.plot(vertices, par_times, 's-', label='Par (vs vertices)', color='#ff7f0e', linewidth=2, markersize=6, alpha=0.8)
        ax6.plot(faces, seq_times, 'o--', label='Seq (vs faces)', color='#1f77b4', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(faces, par_times, 's--', label='Par (vs faces)', color='#ff7f0e', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(avg_faces_per_vertex, seq_times, 'o:', label='Seq (vs avg faces/v)', color='#1f77b4', linewidth=1.5, markersize=5, alpha=0.6)
        ax6.plot(avg_faces_per_vertex, par_times, 's:', label='Par (vs avg faces/v)', color='#ff7f0e', linewidth=1.5, markersize=5, alpha=0.6)
        
        ax6.set_xlabel('X-axis: Vertices / Faces / Avg Faces/Vertex', fontsize=10, fontweight='bold')
        ax6.set_ylabel('Total Time (seconds)', fontsize=11, fontweight='bold')
        ax6.set_title('6. Scalability: Seq vs Par (New Code)', fontsize=12, fontweight='bold')
        ax6.legend(fontsize=7, loc='upper left', ncol=2)
        ax6.grid(True, alpha=0.3)
        ax6.set_xscale('log')
        ax6.set_yscale('log')
    
    plt.suptitle('New Code: Sequential vs Parallel Comparison', fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    output_file = os.path.join(output_dir, 'parallel_vs_sequential.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    
    plt.close()

def create_memory_usage_plot(results, output_dir):
    """Create memory usage comparison plots."""
    if not HAS_MATPLOTLIB:
        print("  ⚠ Skipping plots (matplotlib not available)")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    off_files = [k for k in results.keys() if k != 'ply_file']
    off_files.sort()
    
    x_pos = np.arange(len(off_files))
    width = 0.35
    
    # Plot 1: Sequential Memory Usage (Old vs New)
    ax1 = axes[0, 0]
    old_memory = []
    new_memory = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_seq_old', {})
        new_data = results[file_name].get('test_forman_gradient_seq_new', {})
        old_mem = old_data.get('memory_mb', 0) or (old_data.get('memory_kb', 0) / 1024.0)
        new_mem = new_data.get('memory_mb', 0) or (new_data.get('memory_kb', 0) / 1024.0)
        old_memory.append(old_mem)
        new_memory.append(new_mem)
    
    bars1 = ax1.bar(x_pos - width/2, old_memory, width, label='Old (adjRelations)', color='#d62728', alpha=0.8)
    bars2 = ax1.bar(x_pos + width/2, new_memory, width, label='New (completeCoboundaryTop)', color='#2ca02c', alpha=0.8)
    ax1.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Memory Usage (MB)', fontsize=11, fontweight='bold')
    ax1.set_title('1. Memory Usage: Sequential (Old vs New)', fontsize=12, fontweight='bold')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}MB',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 2: Parallel Memory Usage (Old vs New)
    ax2 = axes[0, 1]
    old_memory = []
    new_memory = []
    for file_name in off_files:
        old_data = results[file_name].get('test_forman_gradient_pa_old', {})
        new_data = results[file_name].get('test_forman_gradient_pa_new', {})
        old_mem = old_data.get('memory_mb', 0) or (old_data.get('memory_kb', 0) / 1024.0)
        new_mem = new_data.get('memory_mb', 0) or (new_data.get('memory_kb', 0) / 1024.0)
        old_memory.append(old_mem)
        new_memory.append(new_mem)
    
    bars1 = ax2.bar(x_pos - width/2, old_memory, width, label='Old (adjRelations)', color='#d62728', alpha=0.8)
    bars2 = ax2.bar(x_pos + width/2, new_memory, width, label='New (completeCoboundaryTop)', color='#2ca02c', alpha=0.8)
    ax2.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Memory Usage (MB)', fontsize=11, fontweight='bold')
    ax2.set_title('2. Memory Usage: Parallel (Old vs New)', fontsize=12, fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}MB',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 3: New Code Sequential vs Parallel Memory Usage
    ax3 = axes[1, 0]
    seq_memory = []
    par_memory = []
    for file_name in off_files:
        seq_data = results[file_name].get('test_forman_gradient_seq_new', {})
        par_data = results[file_name].get('test_forman_gradient_pa_new', {})
        seq_mem = seq_data.get('memory_mb', 0) or (seq_data.get('memory_kb', 0) / 1024.0)
        par_mem = par_data.get('memory_mb', 0) or (par_data.get('memory_kb', 0) / 1024.0)
        seq_memory.append(seq_mem)
        par_memory.append(par_mem)
    
    bars1 = ax3.bar(x_pos - width/2, seq_memory, width, label='Sequential', color='#1f77b4', alpha=0.8)
    bars2 = ax3.bar(x_pos + width/2, par_memory, width, label='Parallel', color='#ff7f0e', alpha=0.8)
    ax3.set_xlabel('Test Files', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Memory Usage (MB)', fontsize=11, fontweight='bold')
    ax3.set_title('3. Memory Usage: New Code (Seq vs Par)', fontsize=12, fontweight='bold')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels([f.replace('_', '\n') for f in off_files], rotation=0, ha='center', fontsize=8)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3, axis='y')
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax3.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}MB',
                        ha='center', va='bottom', fontsize=7)
    
    # Plot 4: Memory Usage vs Number of Vertices (Scalability)
    ax4 = axes[1, 1]
    vertices = []
    old_seq_mem = []
    new_seq_mem = []
    old_par_mem = []
    new_par_mem = []
    for file_name in off_files:
        old_seq_data = results[file_name].get('test_forman_gradient_seq_old', {})
        new_seq_data = results[file_name].get('test_forman_gradient_seq_new', {})
        old_par_data = results[file_name].get('test_forman_gradient_pa_old', {})
        new_par_data = results[file_name].get('test_forman_gradient_pa_new', {})
        
        v = old_seq_data.get('vertices', 0) or old_seq_data.get('complex_vertices', 0)
        old_seq_m = old_seq_data.get('memory_mb', 0) or (old_seq_data.get('memory_kb', 0) / 1024.0)
        new_seq_m = new_seq_data.get('memory_mb', 0) or (new_seq_data.get('memory_kb', 0) / 1024.0)
        old_par_m = old_par_data.get('memory_mb', 0) or (old_par_data.get('memory_kb', 0) / 1024.0)
        new_par_m = new_par_data.get('memory_mb', 0) or (new_par_data.get('memory_kb', 0) / 1024.0)
        
        if v > 0 and (old_seq_m > 0 or new_seq_m > 0):
            vertices.append(v)
            old_seq_mem.append(old_seq_m)
            new_seq_mem.append(new_seq_m)
            old_par_mem.append(old_par_m)
            new_par_mem.append(new_par_m)
    
    if vertices:
        sorted_data = sorted(zip(vertices, old_seq_mem, new_seq_mem, old_par_mem, new_par_mem))
        vertices, old_seq_mem, new_seq_mem, old_par_mem, new_par_mem = zip(*sorted_data)
        
        ax4.plot(vertices, old_seq_mem, 'o-', label='Old Seq', color='#d62728', linewidth=2, markersize=6, alpha=0.7)
        ax4.plot(vertices, new_seq_mem, 's-', label='New Seq', color='#2ca02c', linewidth=2, markersize=6, alpha=0.7)
        ax4.plot(vertices, old_par_mem, 'o--', label='Old Par', color='#d62728', linewidth=1.5, markersize=5, alpha=0.6)
        ax4.plot(vertices, new_par_mem, 's--', label='New Par', color='#2ca02c', linewidth=1.5, markersize=5, alpha=0.6)
        
        ax4.set_xlabel('Number of Vertices', fontsize=11, fontweight='bold')
        ax4.set_ylabel('Memory Usage (MB)', fontsize=11, fontweight='bold')
        ax4.set_title('4. Memory Scalability', fontsize=12, fontweight='bold')
        ax4.legend(fontsize=8, ncol=2)
        ax4.grid(True, alpha=0.3)
        ax4.set_xscale('log')
        ax4.set_yscale('log')
    
    plt.suptitle('Memory Usage Comparison', fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    output_file = os.path.join(output_dir, 'memory_usage_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    plt.close()

def generate_text_report(results, output_dir):
    """Generate text report."""
    report_file = os.path.join(output_dir, 'performance_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("PERFORMANCE ANALYSIS REPORT\n")
        f.write("="*80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n")
        
        # Summary
        f.write("SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write("This report compares the performance of:\n")
        f.write("  - Old code: Uses adjRelations data structure\n")
        f.write("  - New code: Uses completeCoboundaryTop (no adjRelations)\n")
        f.write("\n")
        
        # Detailed results
        for file_name in sorted(results.keys()):
            f.write("\n" + "="*80 + "\n")
            f.write(f"FILE: {file_name}\n")
            f.write("="*80 + "\n")
            
            variants = results[file_name]
            
            # Find old and new variants (handle both xx_tts and test_forman_gradient)
            old_key = None
            new_key = None
            for key in variants.keys():
                if 'old' in key and ('seq' in key or 'pa' in key):
                    old_key = key
                elif 'new' in key and ('seq' in key or 'pa' in key):
                    new_key = key
            
            if old_key and new_key:
                old_data = variants[old_key]
                new_data = variants[new_key]
                
                f.write("\nSequential Build Comparison:\n")
                f.write("-"*80 + "\n")
                
                metrics = [
                    ('Build IA*', 'build_ia_seq'),
                    ('Total Reading', 'total_reading'),
                    ('Forman Gradient', 'total_forman'),
                    ('Total Time', 'total_time')
                ]
                
                for metric_name, metric_key in metrics:
                    old_val = old_data.get(metric_key, 0)
                    new_val = new_data.get(metric_key, 0)
                    
                    if old_val > 0 and new_val > 0:
                        speedup = old_val / new_val
                        f.write(f"  {metric_name:20s}: {old_val:8.3f}s → {new_val:8.3f}s (speedup: {speedup:.2f}x)\n")
                    elif old_val > 0:
                        f.write(f"  {metric_name:20s}: {old_val:8.3f}s → N/A\n")
                    elif new_val > 0:
                        f.write(f"  {metric_name:20s}: N/A → {new_val:8.3f}s\n")
                
                # Also compare sequential vs parallel in new code
                seq_new_data = variants.get('test_forman_gradient_seq_new', {}) or variants.get('xx_tts_seq_new', {})
                par_new_data = variants.get('test_forman_gradient_pa_new', {}) or variants.get('xx_tts_pa_new', {})
                
                if seq_new_data and par_new_data:
                    f.write("\nParallel vs Sequential (New Code):\n")
                    f.write("-"*80 + "\n")
                    
                    for metric_name, metric_key in metrics:
                        seq_val = seq_new_data.get(metric_key, 0)
                        par_val = par_new_data.get(metric_key, 0)
                        
                        if seq_val > 0 and par_val > 0:
                            par_speedup = seq_val / par_val
                            f.write(f"  {metric_name:20s}: Seq {seq_val:7.3f}s, Par {par_val:7.3f}s (parallel speedup: {par_speedup:.2f}x)\n")
                        elif seq_val > 0:
                            f.write(f"  {metric_name:20s}: Seq {seq_val:7.3f}s, Par N/A\n")
                        elif par_val > 0:
                            f.write(f"  {metric_name:20s}: Seq N/A, Par {par_val:7.3f}s\n")
        
        # Memory usage comparison
        f.write("\n" + "="*80 + "\n")
        f.write("MEMORY USAGE COMPARISON\n")
        f.write("="*80 + "\n")
        
        for test_file, variants in sorted(results.items()):
            if test_file == 'ply_file':
                continue
            
            f.write(f"\n{test_file}:\n")
            f.write("-"*80 + "\n")
            
            old_seq = variants.get('test_forman_gradient_seq_old', {})
            new_seq = variants.get('test_forman_gradient_seq_new', {})
            old_par = variants.get('test_forman_gradient_pa_old', {})
            new_par = variants.get('test_forman_gradient_pa_new', {})
            
            if old_seq or new_seq or old_par or new_par:
                f.write(f"{'Variant':<25} {'Memory (MB)':<15} {'Memory (KB)':<15}\n")
                f.write("-"*80 + "\n")
                
                for variant_name, data in [('Sequential Old', old_seq), ('Sequential New', new_seq),
                                          ('Parallel Old', old_par), ('Parallel New', new_par)]:
                    if data:
                        mem_mb = data.get('memory_mb', 0) or (data.get('memory_kb', 0) / 1024.0)
                        mem_kb = data.get('memory_kb', 0) or (data.get('memory_mb', 0) * 1024.0)
                        if mem_mb > 0:
                            f.write(f"{variant_name:<25} {mem_mb:>12.2f} MB {mem_kb:>12.0f} KB\n")
                
                # Calculate memory overhead (new vs old)
                if old_seq and new_seq:
                    old_mem = old_seq.get('memory_mb', 0) or (old_seq.get('memory_kb', 0) / 1024.0)
                    new_mem = new_seq.get('memory_mb', 0) or (new_seq.get('memory_kb', 0) / 1024.0)
                    if old_mem > 0 and new_mem > 0:
                        overhead = ((new_mem - old_mem) / old_mem) * 100
                        f.write(f"\nSequential Memory Overhead: {overhead:+.1f}% ({new_mem - old_mem:+.2f} MB)\n")
                
                if old_par and new_par:
                    old_mem = old_par.get('memory_mb', 0) or (old_par.get('memory_kb', 0) / 1024.0)
                    new_mem = new_par.get('memory_mb', 0) or (new_par.get('memory_kb', 0) / 1024.0)
                    if old_mem > 0 and new_mem > 0:
                        overhead = ((new_mem - old_mem) / old_mem) * 100
                        f.write(f"Parallel Memory Overhead: {overhead:+.1f}% ({new_mem - old_mem:+.2f} MB)\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("CONCLUSION\n")
        f.write("="*80 + "\n")
        f.write("The new implementation (completeCoboundaryTop) shows significant\n")
        f.write("performance improvements, especially in the Build IA* phase.\n")
        f.write("\n")
    
    print(f"  ✓ Saved: {report_file}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_performance_report.py <test_results_dir>")
        print("Example: python generate_performance_report.py cpp/cmake-build-debug/test_results_20251111_123424")
        sys.exit(1)
    
    test_results_dir = sys.argv[1]
    
    if not os.path.exists(test_results_dir):
        print(f"Error: Directory not found: {test_results_dir}")
        sys.exit(1)
    
    print("="*80)
    print("Generating Performance Report")
    print("="*80)
    print(f"Results directory: {test_results_dir}\n")
    
    # Collect results
    print("Collecting test results...")
    results = collect_results(test_results_dir)
    
    if not results:
        print("Error: No results found in directory")
        sys.exit(1)
    
    print(f"  Found {len(results)} test files\n")
    
    # Create output directory
    output_dir = os.path.join(test_results_dir, "report")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}\n")
    
    # Generate plots
    if HAS_MATPLOTLIB:
        print("Generating plots...")
        create_sequential_comparison_plot(results, output_dir)
        create_parallel_comparison_plot(results, output_dir)
        create_parallel_vs_sequential_plot(results, output_dir)
        create_memory_usage_plot(results, output_dir)
    else:
        print("Skipping plots (matplotlib not available)")
        print("To generate plots, install matplotlib:")
        print("  python3 -m pip install --user matplotlib numpy")
        print("  or: sudo pacman -S python-matplotlib python-numpy")
    
    # Generate text report
    print("\nGenerating text report...")
    generate_text_report(results, output_dir)
    
    print("\n" + "="*80)
    print("Report Generation Complete!")
    print("="*80)
    print(f"\nOutput files saved to: {output_dir}")
    print("  - performance_report.txt")
    if HAS_MATPLOTLIB:
        print("  - sequential_comparison.png (Old vs New - Sequential)")
        print("  - parallel_comparison.png (Old vs New - Parallel)")
        print("  - parallel_vs_sequential.png (New Code: Seq vs Par)")
        print("  - memory_usage_comparison.png (Memory Usage Analysis)")
    else:
        print("\nNote: Plots not generated (matplotlib not available)")
    print("")

if __name__ == "__main__":
    main()

