import open3d as o3d
import numpy as np
# Reads point clouds from files pcd1 is bottom linear, pcd2 is bottom diagonal and pcd3 is top linear
pcd1 = o3d.io.read_point_cloud("test_cloud1.ply")
pcd2 = o3d.io.read_point_cloud("test_cloud3.ply")
pcd3 = o3d.io.read_point_cloud("test_cloud2.ply")
print(pcd1, pcd2, pcd3)

# Create frames that show the cameras in the 3D space
camera1 = o3d.geometry.TriangleMesh.create_coordinate_frame()
camera2 = o3d.geometry.TriangleMesh.create_coordinate_frame()
camera3 = o3d.geometry.TriangleMesh.create_coordinate_frame()
origin_frame = o3d.geometry.TriangleMesh.create_coordinate_frame()
origin_frame.scale(0.3, center=origin_frame.get_center())
camera1.scale(0.3, center=camera1.get_center())
camera2.scale(0.3, center=camera2.get_center())
camera3.scale(0.3, center=camera3.get_center())


# Removes any point further than the radius from the camera
def filter_point_cloud(pointcloud, radius):
    points = np.asarray(pointcloud.points)
    center = np.array([0, 0, 0])
    distances = np.linalg.norm(points - center, axis=1)
    pointcloud.points = o3d.utility.Vector3dVector(points[distances <= radius])
    return pointcloud


# Creates the transformation for a 3D object based on angles and translation
def transformation_matrix(geometry, angles, translation):
    T = np.eye(4)
    T[:3, :3] = geometry.get_rotation_matrix_from_xyz(angles)
    T[0, 3] = translation[0]
    T[1, 3] = translation[1]
    T[2, 3] = translation[2]
    return T


filtered_pcd1 = filter_point_cloud(pcd1, 2)
filtered_pcd2 = filter_point_cloud(pcd2, 2)
filtered_pcd3 = filter_point_cloud(pcd3, 2)
extra_filtered_pcd1 = filtered_pcd1
extra_filtered_pcd2 = filtered_pcd2
points = np.asarray(extra_filtered_pcd1.points)
center = np.array([0, 0, 0])
distances = np.linalg.norm(points - center, axis=1)
extra_filtered_pcd1.points = o3d.utility.Vector3dVector(points[distances >= 1.3])
points = np.asarray(extra_filtered_pcd2.points)
center = np.array([0, 0, 0])
distances = np.linalg.norm(points - center, axis=1)
extra_filtered_pcd2.points = o3d.utility.Vector3dVector(points[distances <= 1.3])

o3d.visualization.draw_geometries([extra_filtered_pcd2])
#o3d.io.write_point_cloud("filtered_test1.ply", filtered_pcd1)
#o3d.io.write_point_cloud("filtered_test2.ply", filtered_pcd2)
#o3d.io.write_point_cloud("filtered_test3.ply", filtered_pcd3)
#o3d.io.write_point_cloud("filtered_test4.ply", extra_filtered_pcd1)
o3d.io.write_point_cloud("filtered_test5.ply", extra_filtered_pcd2)