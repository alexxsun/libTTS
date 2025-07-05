'''
test tree detection

this is the runable version Mar-19-2025
use it.
simple method

- denoise by knn distance
- clustering
'''

import os
from collections import defaultdict

import numpy as np
# from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN

from plyfile import PlyData


def dbscan_clustering(pts_np, th_db_eps=0.1, th_db_min_samples=10, set_2d=False):
    # pts: xyzh...
    if set_2d:
        pts_np = pts_np[:, [0, 1]]  # xy
    else:
        pts_np = pts_np[:, [0, 1, 3]]  # xyh
    # print(f"Shape: {pts_np.shape}")
    db = DBSCAN(eps=th_db_eps, min_samples=th_db_min_samples).fit(pts_np)  # 0.1, 10
    labels = db.labels_

    gpids = defaultdict(list)  # label: pids -> pts, to save memory
    # dbscan labels: -1, 0, 1, 2, ..., valid labels start from 0
    for i in range(len(pts_np)):
        l = labels[i]
        if l > -1:
            gpids[l].append(i)
    return labels, gpids


def detect_trees(pts_file, th_low_base=0.0, th_low_length=0.5, th_avg_knn_dist=0.1, out_file=None):
    if not os.path.exists(pts_file):
        print(f"check: {pts_file}")
        return None

    workspace = os.path.dirname(pts_file)

    # read ply by plyfile, get xyzh
    if ".pts" in pts_file:
        pts = np.loadtxt(pts_file)
    else:
        # Read the .ply file
        ply_data = PlyData.read(pts_file)
        # Extract the xyzh coordinates
        vertices = ply_data['vertex']
        pts = np.vstack([vertices['x'], vertices['y'], vertices['z'], vertices['h']]).T

    print(f"pts #: {len(pts)}")
    print(f"pts shape: {pts.shape}")
    # select low height points from xyzh
    pts_hmin = 0.5
    print(f"pts_hmin: {pts_hmin}")
    # find pts between pts_zmin + th_low_base and pts_zmin + th_low_length
    low_pts = pts[(pts[:, 3] > pts_hmin + th_low_base) & (pts[:, 3] < pts_hmin + th_low_base + th_low_length)]  # xyzh
    # low_pts = pts[pts[:, 3] < pts_zmin + th_low_length] # xyzh
    print(f"low pts #: {len(low_pts)}")

    # calculate point wise info of low points
    th_nbr_num = 50
    th_nbr_radius = 1.0

    th_nbr_num = min(th_nbr_num, len(low_pts) - 1)

    # find nearest neighbors for each point in 3D space
    nbrs = NearestNeighbors(n_neighbors=th_nbr_num, algorithm='ball_tree').fit(low_pts[:, [0, 1, 3]])
    distances, indices = nbrs.kneighbors(low_pts[:, [0, 1, 3]])

    # filter points by 
    filtered_low_pts = []
    for i, p in enumerate(low_pts):
        # get neighbors within radius
        nbr_pts = list()
        for nbr, d in zip(indices[i], distances[i]):
            if nbr == i:
                continue
            if d < th_nbr_radius:
                nbr_pts.append(low_pts[nbr])
        if len(nbr_pts) < 3:
            continue

        # find nearest 4 neighbors distance
        nbr_dists = distances[i][1:5]
        avg_knn_dist = np.mean(nbr_dists)
        if avg_knn_dist < th_avg_knn_dist:
            filtered_low_pts.append(p)
    filtered_low_pts = np.array(filtered_low_pts)
    print(f"filtered low pts #: {len(filtered_low_pts)}")

    # clustering
    th_db_eps = 0.02
    th_db_min_samples = 5

    labels, gpids = dbscan_clustering(filtered_low_pts, th_db_eps, th_db_min_samples, set_2d=True)
    labels = np.array(labels)
    # 
    labeled_pts = list()
    for i, p in enumerate(filtered_low_pts):
        x, y, z, h, *_ = p
        l = labels[i]
        # TTS-cpp considers valid labels start from 1
        if l >= 0:
            l += 1
        labeled_pts.append((x, y, z, h, l))
        # labeled_pts.append((x, y, z-th_low_base, h, l)) # make sure z is low to support TTS-cpp. to check.
    labeled_pts = np.array(labeled_pts)

    # save xyzl
    # aka tree location points v1
    if out_file is None:
        out_file = f"{workspace}/test_filtered_labeled_xyzl.pts"
    np.savetxt(out_file, labeled_pts[:, [0, 1, 2, 4]], fmt="%.3f")


import sys

if __name__ == "__main__":
    infile = sys.argv[1]  # .ply file
    if len(sys.argv) > 2:
        outfile = sys.argv[2]  # .pts file

    detect_trees(infile, th_low_base=0.5, th_low_length=0.5, th_avg_knn_dist=0.05, out_file=outfile)
    print("done")
