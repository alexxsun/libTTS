#!/usr/bin/env python3
"""
Compare two PLY segmentation files using Intersection over Union (IoU).

This tool compares segmentation results by computing IoU for each label
and overall. IoU = 1.0 means the files are identical.

Usage:
    python cmp_ply.py file1.ply file2.ply

Exit codes:
    0: Files are identical (IoU = 1.0)
    1: Files are different (IoU < 1.0) or error occurred
"""

import sys
import numpy as np
from plyfile import PlyData


def compute_iou(labels1, labels2):
    """
    Compute Intersection over Union for two label arrays.
    
    Args:
        labels1: numpy array of labels from file 1
        labels2: numpy array of labels from file 2
    
    Returns:
        float: IoU value (0.0 to 1.0)
    """
    # Get unique labels from both files
    all_labels = np.unique(np.concatenate([labels1, labels2]))
    
    if len(all_labels) == 0:
        return 1.0  # Both empty, consider identical
    
    # Compute IoU for each label
    ious = []
    for label in all_labels:
        mask1 = (labels1 == label)
        mask2 = (labels2 == label)
        
        intersection = np.sum(mask1 & mask2)
        union = np.sum(mask1 | mask2)
        
        if union == 0:
            # Both have no points with this label
            ious.append(1.0)
        else:
            iou = intersection / union
            ious.append(iou)
    
    # Overall IoU is the average (or could be weighted by label size)
    overall_iou = np.mean(ious)
    return overall_iou, ious, all_labels


def compare_ply_files(file1_path, file2_path):
    """
    Compare two PLY files with segmentation labels.
    
    Args:
        file1_path: Path to first PLY file
        file2_path: Path to second PLY file
    
    Returns:
        tuple: (overall_iou, per_label_info, success)
    """
    try:
        # Read PLY files
        print(f"Reading {file1_path}...")
        plydata1 = PlyData.read(file1_path)
        
        print(f"Reading {file2_path}...")
        plydata2 = PlyData.read(file2_path)
        
        # Extract vertex data
        vertex1 = plydata1['vertex']
        vertex2 = plydata2['vertex']
        
        # Check if 'label' property exists
        if 'label' not in vertex1.data.dtype.names:
            print(f"Error: 'label' property not found in {file1_path}")
            return None, None, False
        
        if 'label' not in vertex2.data.dtype.names:
            print(f"Error: 'label' property not found in {file2_path}")
            return None, None, False
        
        # Extract labels
        labels1 = vertex1['label']
        labels2 = vertex2['label']
        
        # Check if number of points match
        if len(labels1) != len(labels2):
            print(f"Error: Number of points mismatch!")
            print(f"  {file1_path}: {len(labels1)} points")
            print(f"  {file2_path}: {len(labels2)} points")
            return None, None, False
        
        print(f"Comparing {len(labels1)} points...")
        
        # Compute IoU
        overall_iou, per_label_ious, all_labels = compute_iou(labels1, labels2)
        
        # Prepare per-label information
        per_label_info = []
        for label, iou in zip(all_labels, per_label_ious):
            count1 = np.sum(labels1 == label)
            count2 = np.sum(labels2 == label)
            per_label_info.append({
                'label': label,
                'iou': iou,
                'count_file1': count1,
                'count_file2': count2
            })
        
        return overall_iou, per_label_info, True
        
    except FileNotFoundError as e:
        print(f"Error: File not found: {e}")
        return None, None, False
    except Exception as e:
        print(f"Error reading files: {e}")
        import traceback
        traceback.print_exc()
        return None, None, False


def main():
    if len(sys.argv) != 3:
        print("Usage: python cmp_ply.py <file1.ply> <file2.ply>")
        sys.exit(1)
    
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    
    print("=" * 60)
    print("PLY Segmentation Comparison Tool")
    print("=" * 60)
    print()
    
    overall_iou, per_label_info, success = compare_ply_files(file1_path, file2_path)
    
    if not success:
        print("\nComparison failed!")
        sys.exit(1)
    
    # Print per-label results
    print("\nPer-label IoU:")
    print("-" * 60)
    print(f"{'Label':<10} {'IoU':<15} {'Count (File1)':<15} {'Count (File2)':<15}")
    print("-" * 60)
    
    for info in per_label_info:
        match_status = "✓" if info['iou'] == 1.0 else "✗"
        print(f"{info['label']:<10} {info['iou']:<15.6f} {info['count_file1']:<15} {info['count_file2']:<15} {match_status}")
    
    # Print overall result
    print("-" * 60)
    print(f"\nOverall IoU: {overall_iou:.6f}")
    
    if overall_iou == 1.0:
        print("✓ Files are IDENTICAL")
        sys.exit(0)
    else:
        print("✗ Files are DIFFERENT")
        print(f"  Difference: {1.0 - overall_iou:.6f}")
        sys.exit(1)


if __name__ == "__main__":
    main()

