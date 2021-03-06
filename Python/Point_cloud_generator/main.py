import matplotlib.pyplot as plt
import numpy as np
import open3d as o3d


def read_file(file_name):
    with open(file_name) as file:
        file_text = file.read()
        file_text = file_text.replace("[[", "[")
        file_array = file_text.split("[")[1:]
        for i, element in enumerate(file_array) :
            file_array[i] = element.replace("],", "")
        file_array[-1] = file_array[-1].replace("]]", "")

        for i, distances_text in enumerate(file_array):
            file_array[i] = distances_text.split(",")

        distances = np.empty((len(file_array), len(file_array[0])))

        for i, distances_text in enumerate(file_array):
            for j, distance in enumerate(distances_text):
                distances[i][j] = float(distance)

    distances = np.ndarray.transpose(distances)
    return distances


def get_focus_vector(camera_placement, focus):
    focus_vector = np.array([focus[0] - camera_placement[0], focus[1] - camera_placement[1], focus[2] - camera_placement[2]])
    focus_vector_hat = focus_vector / (focus_vector**2).sum()**0.5
    return focus_vector_hat


def get_focus_axis_1(focus_vector):
    spread_radians = np.pi
    rotational_matrix = [
        [np.cos(spread_radians / 2), -np.sin(spread_radians / 2), 0],
        [np.sin(spread_radians / 2), np.cos(spread_radians / 2), 0],
        [0, 0, 0]
    ]
    vector = np.matmul(rotational_matrix, focus_vector)

    if (vector**2).sum()**0.5 == 0 and focus_vector[2] == 1:
        return np.array([-1, 0, 0])

    elif (vector**2).sum()**0.5 == 0 and focus_vector[2] == -1:
        return np.array([1, 0, 0])

    else:
        return vector / (vector**2).sum()**0.5


def get_focus_axis_2(focus_vector, focus_axis):
    return np.cross(focus_vector, focus_axis)


def get_fov_vectors(vector_1, vector_2, focus_vector, fov):
    fov_radians_side = fov[0] * np.pi / 180 / 2
    fov_radians_top = fov[1] * np.pi / 180 / 2

    leftmost = np.sin(fov_radians_side) + (1 - np.cos(fov_radians_side)) / np.cos(fov_radians_side)
    rightmost = -leftmost
    topmost = np.sin(fov_radians_top) + (1 - np.cos(fov_radians_top)) / np.cos(fov_radians_top)
    bottommost = -topmost

    vector_top_left = leftmost * vector_1 + topmost * vector_2 + focus_vector
    vector_bottom_left = leftmost * vector_1 + bottommost * vector_2 + focus_vector
    vector_bottom_right = rightmost * vector_1 + bottommost * vector_2 + focus_vector
    vector_top_right = rightmost * vector_1 + topmost * vector_2 + focus_vector

    return [vector_top_left, vector_bottom_left, vector_bottom_right, vector_top_right]


def get_point_vectors(vector_1, vector_2, focus_vector, fov, resolution):
    fov_radians_side = fov[0] * np.pi / 180 / 2
    fov_radians_top = fov[1] * np.pi / 180 / 2

    leftmost = np.sin(fov_radians_side) + (1 - np.cos(fov_radians_side))/np.cos(fov_radians_side)
    rightmost = -leftmost
    topmost = np.sin(fov_radians_top) + (1 - np.cos(fov_radians_top))/np.cos(fov_radians_top)
    bottommost = -topmost

    side_vectors_spacing = np.linspace(leftmost, rightmost, resolution[0])
    side_vectors_values = []
    for i in side_vectors_spacing:
        side_vectors_values.append(i * vector_1)

    top_vectors_spacing = np.linspace(topmost, bottommost, resolution[1])
    top_vectors_values = []
    for i in top_vectors_spacing:
        top_vectors_values.append(i * vector_2)

    return side_vectors_values, top_vectors_values


def mapping(value, in_min, in_max, out_min, out_max):
    in_span = in_max - in_min
    out_span = out_max - out_min

    scale = float(value - in_min) / float(in_span)

    return out_min + (scale * out_span)


def get_points(distances, side_vectors, top_vectors, focus_vector):
    point_matrix = []
    grey_scale_value = []
    for i, row in enumerate(distances):
        point_matrix.append([])
        grey_scale_value.append([])
        for j, distance in enumerate(row):
            point_matrix[i].append((side_vectors[j] + top_vectors[i] + focus_vector) * (distance*4))
            grey_scale_value[i].append(mapping(distance, 0, 1, 0, 1))

    return point_matrix, grey_scale_value


def show_point_cloud(cameras, fov_vectors, points, grey_scale, scale_down):
    plt.figure()
    ax = plt.axes(projection='3d')
    plt.xlim(-2, 4)
    plt.ylim(-2, 4)
    ax.set_zlim(-2, 2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    #ax.scatter(X_seen, Y_seen, Z_seen, color="green", marker=".")
    for i, placement in enumerate(cameras):
        ax.scatter(placement[0], placement[1], placement[2], color="black")
        ax.plot([placement[0], placement[0] + fov_vectors[i][0][0]], [placement[1], placement[1] + fov_vectors[i][0][1]],
                [placement[2], placement[2] + fov_vectors[i][0][2]], color="blue")
        ax.plot([placement[0], placement[0] + fov_vectors[i][1][0]], [placement[1], placement[1] + fov_vectors[i][1][1]],
                [placement[2], placement[2] + fov_vectors[i][1][2]], color="blue")
        ax.plot([placement[0], placement[0] + fov_vectors[i][2][0]], [placement[1], placement[1] + fov_vectors[i][2][1]],
                [placement[2], placement[2] + fov_vectors[i][2][2]], color="blue")
        ax.plot([placement[0], placement[0] + fov_vectors[i][3][0]], [placement[1], placement[1] + fov_vectors[i][3][1]],
                [placement[2], placement[2] + fov_vectors[i][3][2]], color="blue")

    for i, row in enumerate(points):
        if i % scale_down == 1 or scale_down == 1 or scale_down == 0:
            for j, point in enumerate(row):
                if j % scale_down == 1 or scale_down == 1 or scale_down == 0:
                    print(f"{round(((i*len(row)+j)/(len(points)*len(row)))*100, 3)}% done")
                    ax.scatter(placement[0] + point[0], placement[1] + point[1], placement[2] + point[2], c="b", marker=".")

    plt.gray()
    ax.view_init(30, 220)
    plt.show()


def main():
    distances = read_file("cam2_dist.txt")
    scale_down = 0
    fov = (39, 25)
    cameras = []
    fov_vectors = []
    cameras.append([0, 0, 0])
    focus_vector = get_focus_vector((0, 0, 0.4), (1, 0.4, 0))
    focus_axis_1 = get_focus_axis_1(focus_vector)
    focus_axis_2 = get_focus_axis_2(focus_vector, focus_axis_1)
    fov_vectors.append(get_fov_vectors(focus_axis_1, focus_axis_2, focus_vector, fov))
    side_vectors, top_vectors = get_point_vectors(focus_axis_1, focus_axis_2, focus_vector, fov, np.shape(distances))
    points, grey_scale = get_points(distances, side_vectors, top_vectors, focus_vector)

    #show_point_cloud(cameras, fov_vectors, points, grey_scale, scale_down)
    ply_point = o3d.data.PLYPointCloud()
    print("Load a ply point cloud, print it, and render it")
    ply_point_cloud = o3d.data.PLYPointCloud()
    pcd = o3d.io.read_point_cloud("1.ply")#ply_point_cloud.path)
    print(len(pcd.points))
    for i, row in enumerate(points[::4]):
        for j, point in enumerate(row):
            pcd.points[i*len(points[0]) + j] = point
    print(np.asarray(pcd.points))
    o3d.visualization.draw_geometries([pcd],
                                      zoom=0.3412,
                                      front=[4, -0.2125, -0.8795],
                                      lookat=[2.6172, 2.0475, 1.532],
                                      up=[-0.0694, -0.9768, 0.2024])

if __name__ == '__main__':
    main()