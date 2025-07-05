
'''
Object based wood leaf separation, or downsample method.
Version: mar-27-2025

given segmentation results (xyzl), e.g., from TTS results.
for each segment, 
    remove the points that are far from the segment, e.g., using avg knn distance.

'''
# try to import the C++ module first, if it exists.
# If the C++ module is not compiled, this will raise an ImportError.
# and set a marker to indicate that the C++ module is not available.
# some functions may not be available if the C++ module is not compiled.

cpp_module_available = False
try:
    from ._libtts import oversegment_tree as _oversegment_tree_cpp

except ImportError:
    # If the import fails, it means the C++ module is not compiled or not found.
    print("Error: C++ module 'libtts' not found. Using Python-only implementation.")

# If the import is successful, set the marker to True.
cpp_module_available = True

import sys
import numpy as np

from plyfile import PlyData, PlyElement

from sklearn.neighbors import NearestNeighbors

import multiprocessing as mp

def select_pts_per_label(pts, th_nbrs = 4, th_dis = 0.1):
        pts = np.array(pts)
        pts = pts[:, :3]
        if len(pts) < th_nbrs:
            return []
        nbrs = NearestNeighbors(n_neighbors=th_nbrs, algorithm='ball_tree').fit(pts)
        distances, indices = nbrs.kneighbors(pts)

        selected_pts = list()
        for i in range(len(pts)):
            tmp_dis = distances[i][1:]
            if np.mean(tmp_dis) < th_dis:
                selected_pts.append(pts[i])
        return selected_pts

def extract_woody_points(infile, th_avg_dis=0.1):
    """
    infile: oversegmentation results, xyzl
    """

    file_ext = infile.split('.')[-1]
    if file_ext not in ['pts', 'ply']:
        print(f"Error: file extension {file_ext} not supported.")
        return 
    if file_ext == 'pts':
        pts_xyzl = np.loadtxt(infile) # xyzl

    if file_ext == 'ply':
        plydata = PlyData.read(infile)
        pts_xyzl = np.array([[x, y, z, l] for x, y, z, l in zip(plydata['vertex']['x'], plydata['vertex']['y'], plydata['vertex']['z'], plydata['l']['pt_label'])])   

    pts_num = len(pts_xyzl)
    print(f"pts #: {pts_num}")

    gpts = dict()# label: pts
    for pid in range(len(pts_xyzl)):
        x, y, z, l = pts_xyzl[pid]
        if l not in gpts:
            gpts[l] = list()
        gpts[l].append([x, y, z, l])
    print(f"labels #: {len(gpts)}")


    num_processes = 8
    
    with mp.Pool(processes=num_processes) as pool:
        # call the function for each label in parallel
        results = pool.starmap(select_pts_per_label, [(pts, 4, th_avg_dis) for l, pts in gpts.items()])
    selected_pts = [pt for pts in results for pt in pts]
    print(f"selected_pts #: {len(selected_pts)}")

    outfl = f"{infile[:-4]}_woody.pts"
    np.savetxt(outfl, selected_pts, fmt="%.3f")
    print(f"Done. Saved to {outfl}")


def downsample_points(infile, th_avg_dis=0.1):
    # use cpp function to process infile
    if not cpp_module_available:
        print("C++ module not available, using Python implementation.")
        extract_woody_points(infile, th_avg_dis)
        return
    print(f"Using C++ implementation to process {infile} with th_avg_dis: {th_avg_dis}")
    # Call the C++ function
    overseg_file = _oversegment_tree_cpp(infile, th_avg_dis)

    # then process the overseg_file
    print(f"Oversegmentation results saved to {overseg_file}")
    extract_woody_points(overseg_file, th_avg_dis)
    

if __name__ == "__main__":
    infile = sys.argv[1]

    th_avg_dis = 0.1
    if len(sys.argv) > 2:
        th_avg_dis = float(sys.argv[2])
    print(f"th_avg_dis: {th_avg_dis}")

    extract_woody_points(infile, th_avg_dis)
    print(f"\nDone")