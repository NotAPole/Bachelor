"""
README

@author: Pål-André Furnes
This project is developed to be used for proof-of-concept in my bachelor thesis.
A pallet is represented by an array of 3D points.
A camera is represented by a single point in space facing a specified direction with a specified field of view.
The resulting 3D-plot shows green and red points. Green means seen, red means unseen.
In fov-check-mode the resulting plot will show points within the camera's fov as green, regardless if it can be seen
or not. This mode is used to verify camera placement.
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing


class Camera:
    def __init__(self, x, y, z, fov, focus, range):
        self.x = x
        self.y = y
        self.z = z
        self.fov = fov
        self.focus = focus
        self.min_range = range[0]
        self.max_range = range[1]

    def check_fov(self, point, possible_obstructions, planar_equations, step, tC, tP):
        x, y, z = self.calculate_coordinates(point, step)
        confirmed_obstructions = self.check_for_obstructions(planar_equations, possible_obstructions, x, y, z, point, step, tC, tP)

        return confirmed_obstructions

    def calculate_coordinates(self, point, step):
        x = np.linspace(self.x, point[0], step)
        y = np.linspace(self.y, point[1], step)
        z = np.linspace(self.z, point[2], step)

        return x, y, z

    @staticmethod
    def check_for_obstructions(equations, surfaces, x, y, z, point, step, threshold_crossing, threshold_point):
        for j in range(step):
            for k, l in enumerate(equations):
                planar_result = l[0] * (x[j] - l[1]) + l[2] * (y[j] - l[3]) + l[4] * (z[j] - l[5])
                crossing_plane = abs(planar_result) < threshold_crossing
                if crossing_plane:
                    corners_x = (np.min([surfaces[k][0][0], surfaces[k][1][0], surfaces[k][2][0], surfaces[k][3][0]]),
                                 np.max([surfaces[k][0][0], surfaces[k][1][0], surfaces[k][2][0], surfaces[k][3][0]]))
                    corners_y = (np.min([surfaces[k][0][1], surfaces[k][1][1], surfaces[k][2][1], surfaces[k][3][1]]),
                                 np.max([surfaces[k][0][1], surfaces[k][1][1], surfaces[k][2][1], surfaces[k][3][1]]))
                    corners_z = (np.min([surfaces[k][0][2], surfaces[k][1][2], surfaces[k][2][2], surfaces[k][3][2]]),
                                 np.max([surfaces[k][0][2], surfaces[k][1][2], surfaces[k][2][2], surfaces[k][3][2]]))
                    confirmed_obstruction = (corners_x[0] - threshold_point) <= x[j] <= (corners_x[1] + threshold_point) and \
                             (corners_y[0] - threshold_point) <= y[j] <= (corners_y[1] + threshold_point) and \
                             (corners_z[0] - threshold_point) <= z[j] <= (corners_z[1] + threshold_point) and \
                             not (point[0] - threshold_point <= x[j] <= point[0] + threshold_point and
                                  point[1] - threshold_point <= y[j] <= point[1] + threshold_point and
                                  point[2] - threshold_point <= z[j] <= point[2] + threshold_point)

                    if confirmed_obstruction:
                        return True

        return False

    def check_points_in_fov(self, points):
        focus_vector = self.get_focus_vector()
        focus_axis_1 = self.get_focus_axis_1(focus_vector)
        focus_axis_2 = get_focus_axis_2(focus_vector, focus_axis_1)
        fov_vectors = get_fov_vectors(focus_axis_1, focus_axis_2, focus_vector, self.fov)
        unit_vectors = []
        for i in range(len(fov_vectors)-1):
            normal_vector = np.cross(fov_vectors[i], fov_vectors[i+1])
            unit_vector = normal_vector / (normal_vector**2).sum()**0.5
            unit_vectors.append(unit_vector)
        normal_vector = np.cross(fov_vectors[-1], fov_vectors[0])
        unit_vector = normal_vector / (normal_vector ** 2).sum() ** 0.5
        unit_vectors.append(unit_vector)

        points_outside_fov = []
        for i in points:
            dist = np.sqrt((i[0] - self.x) ** 2 + (i[1] - self.y)**2 + (i[2] - self.z) ** 2)
            if self.min_range <= dist <= self.max_range:
                for j in unit_vectors:
                    planar_result = np.dot(np.array([i[0], i[1], i[2]]) - np.array([self.x, self.y, self.z]), j)
                    if planar_result > 0:
                        points_outside_fov.append(i)
                        break
                    else:
                        pass
            else:
                points_outside_fov.append(i)

        possible_points = []
        for i in points:
            if i not in points_outside_fov:
                possible_points.append(i)

        return possible_points, points_outside_fov, unit_vectors

    def get_focus_vector(self):
        vector = np.array([self.focus[0]-self.x, self.focus[1]-self.y, self.focus[2]-self.z])
        unit_vector = vector / np.linalg.norm(vector)

        return unit_vector

    @staticmethod
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

    @staticmethod
    def get_spread_vectors_z(focus_vector, rotation_points, spread_degrees): # NOT USED
        four_dimension_focus_vector = np.array([[focus_vector[0]],
                                          [focus_vector[1]],
                                          [focus_vector[2]],
                                          [1]])
        spread_radians = spread_degrees*np.pi/180

        transpose_matrix_1 = np.array([[1, 0, 0, -rotation_points[0][0]],
                                     [0, 1, 0, -rotation_points[0][1]],
                                     [0, 0, 1, -rotation_points[0][2]],
                                     [0, 0, 0, 1]])
        inv_transpose_matrix_1 = np.linalg.inv(transpose_matrix_1)

        rotational_vector = np.array([rotation_points[0][0] - rotation_points[1][0], rotation_points[0][1] - rotation_points[1][1], rotation_points[0][2] - rotation_points[1][2]])
        unit_rotational_vector = rotational_vector / np.linalg.norm(rotational_vector)
        a, b, c = unit_rotational_vector[0], unit_rotational_vector[1], unit_rotational_vector[2]
        d = np.sqrt(b**2 + c**2)
        rotational_matrix_x = np.array([[1, 0, 0, 0],
                                      [0, 1, 0, 0],
                                      [0, 0, 1, 0],
                                      [0, 0, 0, 1]])
        inv_rotational_matrix_x = np.linalg.inv(rotational_matrix_x)
        rotational_matrix_y = np.array([[d, 0, -a, 0],
                                      [0, 1, 0, 0],
                                      [a, 0, d, 0],
                                      [0, 0, 0, 1]])
        inv_rotational_matrix_y = np.linalg.inv(rotational_matrix_y)
        rotational_matrix_z = np.array([[np.cos(spread_radians/2), np.sin(spread_radians/2), 0, 0],
                                      [-np.sin(spread_radians/2), np.cos(spread_radians/2), 0, 0],
                                      [0, 0, 1, 0],
                                      [0, 0, 0, 1]])
        rotated_vector = inv_transpose_matrix_1@inv_rotational_matrix_x@inv_rotational_matrix_y@rotational_matrix_z@rotational_matrix_y@rotational_matrix_x@transpose_matrix_1@four_dimension_focus_vector
        rotated_vector = [rotated_vector[0, 0], rotated_vector[1, 0], rotated_vector[2, 0]]
        unit_rotated_vector = rotated_vector / np.linalg.norm(rotated_vector)

        return unit_rotated_vector


class Pallet:
    def __init__(self, faces):
        self.faces = faces
        self.equations = get_planar_equations(faces)


def show_plots(seen_points, obstructed_points, cameras, unit_vectors=None):
    X_seen, Y_seen, Z_seen = [], [], []
    for j in seen_points:
        X_seen.append(j[0])
        Y_seen.append(j[1])
        Z_seen.append(j[2])

    X_obs, Y_obs, Z_obs = [], [], []
    for j in obstructed_points:
        X_obs.append(j[0])
        Y_obs.append(j[1])
        Z_obs.append(j[2])

    plt.figure()
    ax = plt.axes(projection='3d')
    plt.xlim(-40, 100)
    plt.ylim(-40, 100)
    ax.set_zlim(-50, 50)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.scatter(X_seen, Y_seen, Z_seen, color="green", marker=".")
    ax.scatter(X_obs, Y_obs, Z_obs, color="red", marker=".")
    for i in cameras:
        ax.scatter(i.x, i.y, i.z, color="black")
        fov_vectors = get_camera_fov(i)
        ax.plot([i.x, i.x + fov_vectors[0][0] * 75], [i.y, i.y + fov_vectors[0][1] * 75],
                [i.z, i.z + fov_vectors[0][2] * 75], color="blue")
        ax.plot([i.x, i.x + fov_vectors[1][0] * 75], [i.y, i.y + fov_vectors[1][1] * 75],
                [i.z, i.z + fov_vectors[1][2] * 75], color="blue")
        ax.plot([i.x, i.x + fov_vectors[2][0] * 75], [i.y, i.y + fov_vectors[2][1] * 75],
                [i.z, i.z + fov_vectors[2][2] * 75], color="blue")
        ax.plot([i.x, i.x + fov_vectors[3][0] * 75], [i.y, i.y + fov_vectors[3][1] * 75],
                [i.z, i.z + fov_vectors[3][2] * 75], color="blue")
        if unit_vectors is not None:
            for j in unit_vectors:
                ax.plot([i.x, i.x + j[0] * 10], [i.y, i.y + j[1] * 10],
                        [i.z, i.z + j[2] * 10], color="orange")

    ax.view_init(30, 220)
    plt.show()


def get_planar_equations(surfaces):
    equations = []
    for j in surfaces:
        vector_1 = [j[1][0] - j[0][0], j[1][1] - j[0][1], j[1][2] - j[0][2]]
        vector_2 = [j[2][0] - j[0][0], j[2][1] - j[0][1], j[2][2] - j[0][2]]
        normal = np.cross(vector_1, vector_2)
        normal_hat = normal / (normal ** 2).sum() ** 0.5
        equations.append(get_plane(j[0], normal_hat))

    return equations


def get_camera_fov(camera):
    focus_vector = camera.get_focus_vector()
    focus_axis_1 = camera.get_focus_axis_1(focus_vector)
    focus_axis_2 = get_focus_axis_2(focus_vector, focus_axis_1)
    fov_vectors = get_fov_vectors(focus_axis_1, focus_axis_2, focus_vector, camera.fov)

    return fov_vectors


def get_plane(point, normal):
    return normal[0], point[0], normal[1], point[1], normal[2], point[2]


def get_focus_axis_2(focus_vector, focus_axis):
    return np.cross(focus_vector, focus_axis)


def get_fov_vectors(vector_1, vector_2, focus_vector, fov):
    left_turn_radians = fov[0] * np.pi / 180 / 2
    top_turn_radians = fov[1] * np.pi / 180 / 2

    leftwards = np.sin(left_turn_radians) + (1-np.cos(left_turn_radians))/np.cos(left_turn_radians)
    rightwards = -leftwards
    upwards = np.sin(top_turn_radians) + (1-np.cos(top_turn_radians))/np.cos(top_turn_radians)
    downwards = -upwards

    vector_top_left = leftwards * vector_1 + upwards * vector_2 + focus_vector
    vector_bottom_left = leftwards * vector_1 + downwards * vector_2 + focus_vector
    vector_bottom_right = rightwards * vector_1 + downwards * vector_2 + focus_vector
    vector_top_right = rightwards * vector_1 + upwards * vector_2 + focus_vector

    return [vector_top_left, vector_bottom_left, vector_bottom_right, vector_top_right]


def get_fov_surfaces(point, vectors):
    surfaces = []
    for i in range(len(vectors)-1):
        surfaces.append([point, point + vectors[i], point + vectors[i+1]])
    surfaces.append([point, point + vectors[-1], point + vectors[0]])

    return surfaces


def get_spherical_coordinates(cartesian): # NOT USED
    r = np.sqrt(cartesian[0]**2 + cartesian[1]**2 + cartesian[2]**2)
    theta = np.arctan(cartesian[1]/cartesian[0])
    phi = np.arccos(cartesian[2]/r)

    return [r, theta, phi]


def get_cartesian_coordinates(sphere): # NOT USED
    cartesian_vectors = []
    for i in sphere:
        x = i[0] * np.cos(i[1]) * np.sin(i[2])
        y = i[0] * np.sin(i[1]) * np.sin(i[2])
        z = i[0] * np.cos(i[2])
        cartesian_vectors.append([x, y, z])

    return cartesian_vectors


def get_rotated_vectors(sphere, fov): # NOT USED
    theta_rotation = fov[0]*np.pi/180/2
    phi_rotation = fov[1]*np.pi/180/2
    r, theta, phi = sphere
    top_left_vector = [r, theta+theta_rotation, phi-phi_rotation]
    bottom_left_vector = [r, theta+theta_rotation, phi+phi_rotation]
    bottom_right_vector = [r, theta-theta_rotation, phi+phi_rotation]
    top_right_vector = [r, theta - theta_rotation, phi - phi_rotation]

    return [top_left_vector, bottom_left_vector, bottom_right_vector, top_right_vector]


def check_fov(c, current_camera, pallet, remaining_points, steps, threshold_point, threshold_surface, seen_points=[]):
    camera_progress = 0
    for j in remaining_points:
        pointObstructed = c.check_fov(j, pallet.faces, pallet.equations, steps, threshold_surface, threshold_point)
        if pointObstructed:
            pass
        else:
            seen_points.append(j)
        camera_progress += 1
        print("Camera " + str(current_camera) + " " + str(
            round(camera_progress / len(remaining_points) * 100, 2)) + "% done")

    return seen_points


def read_file(file, additional_points):
    with open(file) as palletFile:
        pallet_surface_text = palletFile.read()
        pallet_surface_text = pallet_surface_text.replace("[", "")
        pallet_surface_text = pallet_surface_text.replace("]", "")
        pallet_surface_text = pallet_surface_text.replace(" ", "")
        pallet_surface_text = pallet_surface_text.split(",")
        pallet_surface_num = []
        initial_points = []
        for i in range(int(len(pallet_surface_text) / 12.0)):
            pallet_surface_temp = []
            for j in range(4):
                initial_points.append([float(pallet_surface_text[i * 12 + j * 3]),
                                        float(pallet_surface_text[i * 12 + 1 + j * 3]),
                                        float(pallet_surface_text[i * 12 + 2 + j * 3])])
                pallet_surface_temp.append([float(pallet_surface_text[i * 12 + j * 3]),
                                          float(pallet_surface_text[i * 12 + 1 + j * 3]),
                                          float(pallet_surface_text[i * 12 + 2 + j * 3])])
            pallet_surface_num.append(pallet_surface_temp)
        for j in initial_points:
            if initial_points.count(j) > 1:
                initial_points.remove(j)

        for j in pallet_surface_num:
            if j[0][0] != j[1][0]:
                for k in np.linspace(j[0][0], j[1][0], additional_points, endpoint=False):
                    initial_points.append([k, j[0][1], j[0][2]])
                for k in np.linspace(j[2][0], j[3][0], additional_points, endpoint=False):
                    initial_points.append([k, j[2][1], j[2][2]])

                if j[1][1] != j[2][1]:
                    for k in np.linspace(j[1][1], j[2][1], additional_points, endpoint=False):
                        initial_points.append([j[1][0], k, j[1][2]])
                    for k in np.linspace(j[3][1], j[0][1], additional_points, endpoint=False):
                        initial_points.append([j[3][0], k, j[3][2]])

                if j[1][2] != j[2][2]:
                    for k in np.linspace(j[1][2], j[2][2], additional_points, endpoint=False):
                        initial_points.append([j[1][0], j[1][1], k])
                    for k in np.linspace(j[3][2], j[0][2], additional_points, endpoint=False):
                        initial_points.append([j[3][0], j[3][1], k])

            if j[0][1] != j[1][1]:
                for k in np.linspace(j[0][1], j[1][1], additional_points, endpoint=False):
                    initial_points.append([j[0][0], k, j[0][2]])
                for k in np.linspace(j[2][1], j[3][1], additional_points, endpoint=False):
                    initial_points.append([j[3][0], k, j[2][2]])

                if j[1][0] != j[2][0]:
                    for k in np.linspace(j[1][0], j[2][0], additional_points, endpoint=False):
                        initial_points.append([k, j[1][1], j[1][2]])
                    for k in np.linspace(j[3][0], j[0][0], additional_points, endpoint=False):
                        initial_points.append([k, j[3][1], j[3][2]])

                if j[1][2] != j[2][2]:
                    for k in np.linspace(j[1][2], j[2][2], additional_points, endpoint=False):
                        initial_points.append([j[1][0], j[1][1], k])
                    for k in np.linspace(j[3][2], j[0][2], additional_points, endpoint=False):
                        initial_points.append([j[3][0], j[3][1], k])

            if j[0][2] != j[1][2]:
                for k in np.linspace(j[0][2], j[1][2], additional_points, endpoint=False):
                    initial_points.append([j[0][0], j[0][1], k])
                for k in np.linspace(j[2][0], j[3][0], additional_points, endpoint=False):
                    initial_points.append([[j[2][0]], j[2][1], k])

                if j[1][0] != j[2][0]:
                    for k in np.linspace(j[1][0], j[2][0], additional_points, endpoint=False):
                        initial_points.append([k, j[1][1], j[1][2]])
                    for k in np.linspace(j[3][0], j[0][0], additional_points, endpoint=False):
                        initial_points.append([k, j[3][1], j[3][2]])

                if j[1][1] != j[2][1]:
                    for k in np.linspace(j[1][1], j[2][1], additional_points, endpoint=False):
                        initial_points.append([j[1][0], k, j[1][2]])
                    for k in np.linspace(j[3][1], j[0][1], additional_points, endpoint=False):
                        initial_points.append([j[3][0], k, j[3][2]])

        for j in initial_points:
            if initial_points.count(j) > 1:
                initial_points.remove(j)

    return pallet_surface_num, initial_points


def main():
    debug = False
    fov_check = True
    use_multi_processing = False
    steps = 1
    additional_points = 20
    threshold_surface = 1/steps*20
    threshold_point = 0.02
    start_time = time.time()
    pallet_surface_num, initial_points = read_file("euro.txt", additional_points)
    if debug:
        cameras = []
        cameras.append(Camera(-10, -10, 0, (50, 36), (0, 0, 0), (0-150)))
        p1 = Pallet(pallet_surface_num)

        points_in_fov, not_possible_points, unit_vectors = cameras[0].check_points_in_fov(initial_points)
        seen_points_camera = check_fov(cameras[0], 1, p1, points_in_fov, steps, threshold_point, threshold_surface)

        obstructed_points = []
        for i in initial_points:
            if i not in seen_points_camera:
                obstructed_points.append(i)
        show_plots(points_in_fov, not_possible_points, cameras, unit_vectors)

    elif fov_check:
        cameras = []
        cameras.append(Camera(-70, -60, 7, (35, 29), (0, 5, 7), (170, 200)))

        points_in_fov, not_possible_points, unit_vectors = cameras[0].check_points_in_fov(initial_points)
        show_plots(points_in_fov, not_possible_points, cameras, unit_vectors)

    else:
        cameras = []
        cameras.append(Camera(-74.7, -125.3, -86, (50, 25), (0, 0, 10), (80, 260)))
        #cameras.append(Camera(116.4, 257.2, -86, (70, 55), (0, 0, 0), (80, 260)))
        #cameras.append(Camera(60, 260, 94.4, (70, 55), (40, 60, 10), (80, 260)))
        #cameras.append(Camera(0, -50, 0, (70, 55), (0, 0, 20), (80, 260)))
        p1 = Pallet(pallet_surface_num)
        if use_multi_processing:
            manager = multiprocessing.Manager()
            seen_points = manager.list()
            process_pool = []
            for i, c in enumerate(cameras):
                remaining_points, obstructed_points = c.check_points_in_fov(initial_points)
                process_pool.append(multiprocessing.Process(target=check_fov, args=(c, i+1, p1, remaining_points, steps, threshold_point, threshold_surface, seen_points)))
                process_pool[i].start()
            for process in process_pool:
                process.join()

            for i in seen_points:
                if seen_points.count(i) > 1:
                    seen_points.remove(i)
        else:
            seen_points = []
            remaining_points = initial_points
            for i, c in enumerate(cameras):
                points_in_fov, _, _ = c.check_points_in_fov(remaining_points)
                seen_points_camera = check_fov(c, i+1, p1, points_in_fov, steps, threshold_point, threshold_surface)
                remaining_points = []
                for j in seen_points_camera:
                    seen_points.append(j)
                for j in initial_points:
                    if j not in seen_points:
                        remaining_points.append(j)

        obstructed_points = []
        for i in initial_points:
            if i not in seen_points:
                obstructed_points.append(i)

        print("Total time used:", time.time() - start_time)

        show_plots(seen_points, obstructed_points, cameras)


if __name__ == "__main__":
    main()
