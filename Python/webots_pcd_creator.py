import numpy as np
import open3d as o3d


def get_calibration_matrix():
    file_data = np.loadtxt('zivid_test/k_matrix.txt').reshape(3, 3)
    return file_data


# Uses the depth array and k matrix to create a point cloud,
# the axis are rotated to correspond to the orientation in webots
def get_point_cloud(depth_image, k_matrix):
    inv_fx = 1.0 / k_matrix[0, 0]
    inv_fy = 1.0 / k_matrix[1, 1]
    ox = k_matrix[0, 2]
    oy = k_matrix[1, 2]
    image_height, image_width = depth_image.shape
    points = np.zeros((image_width * image_height, 3), dtype=np.float32)
    counter = 0
    for y in range(image_height):
        for x in range(image_width):
            dist = depth_image[y, x]
            points[counter, 2] = -np.float32((x - ox) * dist * inv_fx)
            points[counter, 1] = -np.float32((y - oy) * dist * inv_fy)
            points[counter, 0] = np.float32(dist)
            counter += 1
    return points[:counter].astype(np.float32)


# Reads array written from webots, needs the file name and shape of array
def read_array(file_name, x_px, y_px, z):
    if z == 0:
        file_data = np.loadtxt(file_name).reshape(x_px, y_px)
    else:
        file_data = np.loadtxt(file_name).reshape(x_px, y_px, z)
    return file_data


def save_point_cloud(name, pcd):
    o3d.io.write_point_cloud(name, pcd)


# Creates coordinate frames for the origin and cameras and scales them down, so they don't dominate the view
def create_cameras():
    camera1 = o3d.geometry.TriangleMesh.create_coordinate_frame()
    camera2 = o3d.geometry.TriangleMesh.create_coordinate_frame()
    camera3 = o3d.geometry.TriangleMesh.create_coordinate_frame()
    origin_frame = o3d.geometry.TriangleMesh.create_coordinate_frame()
    origin_frame.scale(0.3, center=origin_frame.get_center())
    camera1.scale(0.3, center=camera1.get_center())
    camera2.scale(0.3, center=camera2.get_center())
    camera3.scale(0.3, center=camera3.get_center())
    return origin_frame, camera1, camera2, camera3


def create_transformation(geometry, quat, translation):
    T = np.eye(4)
    T[:3, :3] = geometry.get_rotation_matrix_from_quaternion(quat)
    T[0, 3] = translation[0]
    T[1, 3] = translation[1]
    T[2, 3] = translation[2]
    return T


def filter_point_cloud(pointcloud, radius):
    points = np.asarray(pointcloud.points)
    colours = np.asarray(pointcloud.colors)
    center = np.array([0, 0, 0])
    distances = np.linalg.norm(points - center, axis=1)
    pointcloud.points = o3d.utility.Vector3dVector(points[distances <= radius])
    pointcloud.colors = o3d.utility.Vector3dVector(colours[distances <= radius])
    return pointcloud


# camera resolution
width = 1920
height = 1200

# collects the images from text files

k_matrix = get_calibration_matrix()
dist_image_2 = read_array('zivid_test/cam2_dist.txt', width, height, 0)
point_cloud_2 = get_point_cloud(dist_image_2, k_matrix)
colour_img_2 = read_array('zivid_test/cam2_rgb.txt', width, height, 3)
colour_2 = colour_img_2.reshape(height*width, 3)
dist_image_1 = read_array('zivid_test/cam1_dist.txt', width, height, 0)
point_cloud_1 = get_point_cloud(dist_image_1, k_matrix)
colour_img_1 = read_array('zivid_test/cam1_rgb.txt', width, height, 3)
colour_1 = colour_img_1.reshape(height*width, 3)
dist_image_3 = read_array('zivid_test/cam3_dist.txt', width, height, 0)
point_cloud_3 = get_point_cloud(dist_image_3, k_matrix)
colour_img_3 = read_array('zivid_test/cam3_rgb.txt', width, height, 3)
colour_3 = colour_img_3.reshape(height*width, 3)

# creates the geometries and displays them
pcd1 = o3d.geometry.PointCloud()
pcd2 = o3d.geometry.PointCloud()
pcd3 = o3d.geometry.PointCloud()
pcd1.points = o3d.utility.Vector3dVector(point_cloud_1)
pcd1.colors = o3d.utility.Vector3dVector(colour_1/255)
pcd2.points = o3d.utility.Vector3dVector(point_cloud_2)
pcd2.colors = o3d.utility.Vector3dVector(colour_2/255)
pcd3.points = o3d.utility.Vector3dVector(point_cloud_3)
pcd3.colors = o3d.utility.Vector3dVector(colour_3/255)

origin, cam1, cam2, cam3 = create_cameras()

pcd1 = filter_point_cloud(pcd1, 2.999)
pcd2 = filter_point_cloud(pcd2, 2.999)
pcd3 = filter_point_cloud(pcd3, 2.999)
save_point_cloud('zivid_top_view.ply', pcd1)
save_point_cloud('zivid_bot_view_1.ply', pcd2)
save_point_cloud('zivid_bot_view_2.ply', pcd3)

T1 = create_transformation(pcd1, [0.926121, -5.6e-07, 0.377227, -9.1717e-07], [-0.283, 1.45, 2.81])
cam1.transform(T1)
pcd1.transform(T1)
T2 = create_transformation(pcd2, [0.933011, -0.12941, -0.224145, -0.250005], [-0.455, 2.08, 0.677])
cam2.transform(T2)
pcd2.transform(T2)
T3 = create_transformation(pcd3, [0.254406, 0.25441, -0.033494, 0.932433], [2.25, 0.777, 0.697])
pcd3.transform(T3)
cam3.transform(T3)

o3d.visualization.draw_geometries([pcd1, pcd2, pcd3, cam1, cam2, cam3, origin])
save_point_cloud('zivid_top_view.ply', pcd1)
save_point_cloud('zivid_bot_view_1.ply', pcd2)
save_point_cloud('zivid_bot_view_2.ply', pcd3)


combined_pcd = o3d.geometry.PointCloud()
combined_pcd = pcd1 + pcd2 + pcd3

o3d.io.write_point_cloud("Zivid_combined.ply", combined_pcd)
