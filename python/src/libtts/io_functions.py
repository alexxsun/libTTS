from plyfile import PlyData
import numpy as np

# show summary of ply file: pts number, vertex element names, type etc
def show_ply_point_summary(filepath):
    """
    Reads a PLY file and prints a summary of its point data.

    Args:
        filepath (str): The path to the PLY file.
    """
    try:
        plydata = PlyData.read(filepath)
        
        # show element names and properties
        print(f"Reading PLY file: {filepath}")
        print("Elements in the PLY file:")
        for element in plydata.elements:
            print(f"  - {element.name}: {len(element.data)} points")
            print(f"    Properties:")
            for prop in element.data.dtype.names:
                print(f"      - {prop}: {element.data.dtype[prop]}")
        print("\nPoint data summary:")

        element_names = [element.name for element in plydata.elements]
        if 'vertex' in element_names:
            vertex_data = plydata['vertex'].data
            print(f"Number of points: {len(vertex_data)}")

            # # Print information about the properties (e.g., x, y, z, color)
            # print("Point properties and their data types:")
            # for prop_name in vertex_data.dtype.names:
            #     print(f"  - {prop_name}: {vertex_data.dtype[prop_name]}")

            # Print min/max values for specific properties 
            # if # of pts is relatively small, otherwise it may take too long
            if len(vertex_data) < 1000000 and \
               'x' in vertex_data.dtype.names and \
               'y' in vertex_data.dtype.names and \
               'z' in vertex_data.dtype.names:
                print("\nCoordinate ranges:")
                print(f"  X-range: [{vertex_data['x'].min():.4f}, {vertex_data['x'].max():.4f}]")
                print(f"  Y-range: [{vertex_data['y'].min():.4f}, {vertex_data['y'].max():.4f}]")
                print(f"  Z-range: [{vertex_data['z'].min():.4f}, {vertex_data['z'].max():.4f}]")
        else:
            print(f"No 'vertex' element found in '{filepath}'.")
            print("Available elements:", [el.name for el in plydata.elements])

    except FileNotFoundError:
        print(f"Error: File not found at '{filepath}'")
    except Exception as e:
        print(f"An error occurred: {e}")

def read_points_from_ply(file_path):
    """
    Reads point cloud data from a PLY file.

    Args:
        file_path (str): Path to the PLY file.

    Returns:
        numpy.ndarray: Array of points read from the PLY file.
    """
    ply_data = PlyData.read(file_path)
    points = ply_data['vertex'].data
    return points

def save_points_xyzh_to_ply(points, filename):
    from plyfile import PlyData, PlyElement
    vertex = np.array([(x, y, z, h) for x, y, z, h in points], 
                      dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('h', 'f4')])
    el = PlyElement.describe(vertex, 'vertex')
    PlyData([el]).write(filename)
    print(f"Saved {len(points)} points to {filename}")
    return

def read_points_from_txt(file_path):
    """
    Reads point cloud data from a text file.

    Args:
        file_path (str): Path to the text file.

    Returns:
        numpy.ndarray: Array of points read from the text file.
    """
    import numpy as np
    points = np.loadtxt(file_path, delimiter=',')
    return points