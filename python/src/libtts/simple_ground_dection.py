import cv2
import numpy as np
from scipy.interpolate import griddata
from plyfile import PlyData, PlyElement

import matplotlib.pyplot as plt


def plot_ground_grid(grid_x, grid_y, grid_z, title):
    # Plotting the 2D heatmap
    plt.figure(figsize=(10, 8))
    plt.imshow(grid_z, extent=(grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()), origin='lower', cmap='terrain')
    plt.colorbar(label='Z value')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title(title)
    plt.show()
    plt.close()


# get ground grid
def get_ground_grid(pts_xyz, grid_size=0.1, show_plot=False):
    # divide points into grid, find the lowest point in each grid, smooth the ground as a plane

    pts_xyz = np.array(pts_xyz)  # ensure numpy array
    pts_xyz = pts_xyz[:, :3]  # (n, 3)

    # find the bounding box of the points
    min_x = np.min(pts_xyz[:, 0])
    max_x = np.max(pts_xyz[:, 0])
    min_y = np.min(pts_xyz[:, 1])
    max_y = np.max(pts_xyz[:, 1])
    min_z = np.min(pts_xyz[:, 2])
    max_z = np.max(pts_xyz[:, 2])
    if show_plot:
        print(f"{min_x=}, {max_x=}, {min_y=}, {max_y=}, {min_z=}, {max_z=}")

    # create a grid
    grid_x = np.arange(min_x, max_x, grid_size)
    grid_y = np.arange(min_y, max_y, grid_size)
    grid_x, grid_y = np.meshgrid(grid_x, grid_y)
    grid_z = np.full_like(grid_x, np.nan)  # default to nan
    # grid_z_ptsnum = np.zeros_like(grid_x) # number of points in each grid

    # update the grid_z by checking each point located in which grid and 
    for p in pts_xyz:
        x, y, z = p
        i = int((x - min_x) / grid_size)
        j = int((y - min_y) / grid_size)
        # grid_z_ptsnum[j, i] += 1
        if np.isnan(grid_z[j, i]) or z < grid_z[j, i]:
            grid_z[j, i] = z
    if show_plot:
        plot_ground_grid(grid_x, grid_y, grid_z, "ground grid: initial")

    # remove the outliers (e.g., 1 sigma), especially at the edge
    grid_nonan = grid_z[~np.isnan(grid_z)]
    mean_z = np.mean(grid_nonan)
    std_z = np.std(grid_nonan)
    print(f"mean_z: {mean_z}, std_z: {std_z}")
    grid_z[grid_z > mean_z + 1 * std_z] = np.nan
    grid_z[grid_z < mean_z - 1 * std_z] = np.nan
    if show_plot:
        plot_ground_grid(grid_x, grid_y, grid_z, "ground grid: remove outliers")

    # fill the nan values: 
    nan_mask = np.isnan(grid_z)
    if len(np.where(nan_mask)[0]) > 0:
        print(f"nan #: {len(np.where(nan_mask)[0])} / {grid_x.shape[0] * grid_x.shape[1]}")
        # Get the coordinates of the valid (non-NaN) points
        valid_coords = np.array(np.nonzero(~nan_mask)).T
        valid_values = grid_z[~nan_mask]
        # Get the coordinates of the NaN points
        nan_coords = np.array(np.nonzero(nan_mask)).T
        # Interpolate the NaN values using griddata
        grid_z[nan_mask] = griddata(valid_coords, valid_values, nan_coords, method='linear')  # linear
        # If there are still NaN values (e.g., if they were outside the convex hull of the valid points), you can use nearest-neighbor interpolation as a fallback
        if np.isnan(grid_z).any():
            # grid_z[nan_mask] = griddata(valid_coords, valid_values, nan_coords, method='nearest')
            nan_num = len(np.where(np.isnan(grid_z))[0])
            print(f"outside nan #: {nan_num} / {grid_x.shape[0] * grid_x.shape[1]}")

    if show_plot:
        plot_ground_grid(grid_x, grid_y, grid_z, "ground grid: after fill nan")

    # smooth the ground
    grid_z = cv2.GaussianBlur(grid_z, (5, 5), 0)
    if show_plot:
        plot_ground_grid(grid_x, grid_y, grid_z, "ground grid: after smooth")

    return grid_x, grid_y, grid_z


# Interpolation, using LinearNDInterpolator
from scipy.interpolate import LinearNDInterpolator


def normalize_pts(pts_xyz, grid_width, outfile, show_plot=False):
    grid_x, grid_y, grid_z = get_ground_grid(pts_xyz, grid_size=grid_width, show_plot=show_plot)

    # get ground grid (x,y,z)
    grid_xyz = []
    for i in range(grid_x.shape[0]):
        for j in range(grid_x.shape[1]):
            grid_xyz.append([grid_x[i, j], grid_y[i, j], grid_z[i, j]])
    grid_xyz = np.array(grid_xyz)
    if ".pts" in outfile:
        grid_outfile = outfile.replace(".pts", "_grid.pts")
    elif ".ply" in outfile:
        grid_outfile = outfile.replace(".ply", "_grid.pts")
    else:
        print("output file should be .pts or .ply")
        sys.exit(0)
    # grid pts shouldn't be too large, .pts is enough and easy to read
    np.savetxt(grid_outfile, grid_xyz)
    print(f"ground grid saved to {grid_outfile}")

    # get ground height from grid
    # grid_x, grid_y are meshgrid
    # grid_z is the height
    # all are numpy 2D array
    interp = LinearNDInterpolator(grid_xyz[:, 0:2], grid_xyz[:, 2])

    xs = np.array([p[0] for p in pts_xyz])
    ys = np.array([p[1] for p in pts_xyz])
    hs = interp(xs, ys)

    # xyzh
    out_xyzh = []
    for i in range(len(pts_xyz)):
        x = pts_xyz[i][0]
        y = pts_xyz[i][1]
        z = pts_xyz[i][2]
        h = hs[i]
        out_xyzh.append([x, y, z, z - h])
    out_xyzh = np.array(out_xyzh)

    # write normalization results
    if ".pts" in outfile:
        np.savetxt(outfile, out_xyzh)
    elif ".ply" in outfile:
        # write ply file
        vertex = np.array([tuple(p) for p in out_xyzh], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('h', 'f4')])
        el = PlyElement.describe(vertex, 'vertex')
        PlyData([el]).write(outfile)  # ouput binary ply file
    else:
        print("output file should be .pts or .ply")
        sys.exit(0)
    print(f"normalized points saved to {outfile}")
    return out_xyzh


def classify_pts(pts_xyzh, infile, th_h=0.5):
    # veg pts:
    veg_pts = [p for p in pts_xyzh if p[3] > th_h]  # xyz
    veg_pts = np.array(veg_pts)  # (n, 4)
    if ".pts" in infile:
        # xyzh
        veg_ptsfile = infile.replace(".pts", "_veg.pts")
        np.savetxt(veg_ptsfile, veg_pts)
        print(f"veg points saved to {veg_ptsfile}")
        # xyh
        veg_xyh_ptsfile = infile.replace(".pts", "_veg_xyh.pts")
        np.savetxt(veg_xyh_ptsfile, veg_pts[:, [0, 1, 3]])
        print(f"veg points saved to {veg_xyh_ptsfile}")
        # xyz: for segmentation task
        veg_xyz_ptsfile = infile.replace(".pts", "_veg_xyz.pts")
        np.savetxt(veg_xyz_ptsfile, veg_pts[:, :3])
        print(f"veg points saved to {veg_xyz_ptsfile}")
    elif ".ply" in infile:
        # xyzh
        veg_ptsfile = infile.replace(".ply", "_veg.ply")
        vertex = np.array([tuple(p) for p in veg_pts], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('h', 'f4')])
        el = PlyElement.describe(vertex, 'vertex')
        PlyData([el]).write(veg_ptsfile)
        print(f"veg points saved to {veg_ptsfile}")
        # xyh
        veg_xyh_ptsfile = infile.replace(".ply", "_veg_xyh.ply")
        vertex = np.array([tuple(p) for p in veg_pts[:, [0, 1, 3]]], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4')])
        el = PlyElement.describe(vertex, 'vertex')
        PlyData([el]).write(veg_xyh_ptsfile)
        print(f"veg points saved to {veg_xyh_ptsfile}")
        # xyz: for segmentation task
        veg_xyz_ptsfile = infile.replace(".ply", "_veg_xyz.ply")
        vertex = np.array([tuple(p) for p in veg_pts[:, :3]], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4')])
        el = PlyElement.describe(vertex, 'vertex')
        PlyData([el]).write(veg_xyz_ptsfile)
        print(f"veg points saved to {veg_xyz_ptsfile}")
    else:
        print("input file should be .pts or .ply")
        sys.exit(0)

    # ground pts
    gd_pts = [p for p in pts_xyzh if p[3] <= th_h]  # xyzh
    gd_pts = np.array(gd_pts)
    if ".pts" in infile:
        # xyzh
        gd_ptsfile = infile.replace(".pts", "_gd.pts")
        np.savetxt(gd_ptsfile, gd_pts)
        print(f"ground points saved to {gd_ptsfile}")
        # xyh
        gd_xyh_ptsfile = infile.replace(".pts", "_gd_xyh.pts")
        np.savetxt(gd_xyh_ptsfile, gd_pts[:, [0, 1, 3]])
        print(f"ground points saved to {gd_xyh_ptsfile}")
    elif ".ply" in infile:
        # xyzh
        gd_ptsfile = infile.replace(".ply", "_gd.ply")
        vertex = np.array([tuple(p) for p in gd_pts], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('h', 'f4')])
        el = PlyElement.describe(vertex, 'vertex')
        PlyData([el]).write(gd_ptsfile)
        print(f"ground points saved to {gd_ptsfile}")
        # xyh
        gd_xyh_ptsfile = infile.replace(".ply", "_gd_xyh.ply")
        vertex = np.array([tuple(p) for p in gd_pts[:, [0, 1, 3]]], dtype=[('x', 'f4'), ('y', 'f4'), ('h', 'f4')])
        el = PlyElement.describe(vertex, 'vertex')
        PlyData([el]).write(gd_xyh_ptsfile)
        print(f"ground points saved to {gd_xyh_ptsfile}")
    else:
        print("input file should be .pts or .ply")
        sys.exit(0)
    return


def normalize_pts_with_gd_grid(pts_xyz, gd_grid, outfile):
    # pts_xyz: (n, 3)
    # gd_grid: (m, 3)

    interp = LinearNDInterpolator(gd_grid[:, 0:2], gd_grid[:, 2])

    xs = np.array([p[0] for p in pts_xyz])
    ys = np.array([p[1] for p in pts_xyz])
    hs = interp(xs, ys)

    # xyzh
    out_xyzh = []
    for i in range(len(pts_xyz)):
        x = pts_xyz[i][0]
        y = pts_xyz[i][1]
        z = pts_xyz[i][2]
        h = hs[i]
        out_xyzh.append([x, y, z, z - h])
    out_xyzh = np.array(out_xyzh)

    # write normalization results
    if ".pts" in outfile:
        np.savetxt(outfile, out_xyzh)
    elif ".ply" in outfile:
        # write ply file
        vertex = np.array([tuple(p) for p in out_xyzh], dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('h', 'f4')])
        el = PlyElement.describe(vertex, 'vertex')
        PlyData([el]).write(outfile)  # ouput binary ply file
    else:
        print("output file should be .pts or .ply")
        sys.exit(0)
    print(f"normalized points saved to {outfile}")
    return out_xyzh


import sys

if __name__ == "__main__":

    if sys.argv[-1] == "-v2":
        print(f"normalize pts with provided ground grid\n")
        infile = sys.argv[1]  # .pts file
        gdfile = sys.argv[2]  # .pts file

        pts_xyz = np.loadtxt(infile)

        if ".pts" in gdfile:
            gd_xyz = np.loadtxt(gdfile)
        elif ".ply" in gdfile:
            plydata = PlyData.read(gdfile)
            gd_xyz = np.vstack([plydata['vertex']['x'], plydata['vertex']['y'], plydata['vertex']['z']]).T
        else:
            print("input file should be .pts or .ply")
            sys.exit(0)

        outfile = infile.replace(".pts", "_normalized.pts")
        print(f"output file: {outfile}")

        normalize_pts_with_gd_grid(pts_xyz, gd_xyz, outfile)

        print(f"Done")
        exit(0)

    # todo: add argparse, update/merge the code
    infile = sys.argv[1]  # .pts file
    gd_width = float(sys.argv[2])  # grid width, 0.1m
    th_veg_h = 0.5  # threshold for classifying ground and vegetation points
    if len(sys.argv) > 3:
        th_veg_h = float(sys.argv[3])

    if ".pts" in infile:
        print(f"loading {infile}")
        inpts = np.loadtxt(infile)
        outfile = infile.replace(".pts", "_normalized.pts")
        print(f"output file: {outfile}")
    elif ".ply" in infile:
        print(f"loading {infile}")
        plydata = PlyData.read(infile)
        inpts = np.vstack([plydata['vertex']['x'], plydata['vertex']['y'], plydata['vertex']['z']]).T
        outfile = infile.replace(".ply", "_normalized.ply")
        print(f"output file: {outfile}")
    else:
        print("input file should be .pts or .ply")
        sys.exit(0)

    pts_xyzh = normalize_pts(inpts, gd_width, outfile, show_plot=True)

    print("classifying ground and vegetation points")
    print(f"threshold: {th_veg_h}")
    classify_pts(pts_xyzh, infile, th_veg_h)

    print("done")
    sys.exit(0)
