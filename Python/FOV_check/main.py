"""
README

@author: Pål-André Furnes
This project is developed to be used for proof-of-concept in my bachelor's thesis.
A pallet is represented by an array of 3D points.
A camera is represented by a single point in space facing a specified direction along with a specified field of view.
The resulting 3D-plot shows green and red points. Green means seen, red means unseen.
In fov-check-mode the resulting plot will show points within the camera's fov as blue, regardless if it can actually
be seen or not. This mode is used to verify camera placements where the objective is to see as much of the points as
possible.
"""

import time
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing


class Queue:
    def __init__(self):
        self.queue = []

    # Check to see whether or not queue is empty
    def isEmpty(self) -> bool:
        return True if len(self.queue) == 0 else False

    # Return the first element of the queue
    def front(self):
        return self.queue[0]

    # Return the last element of the queue
    def rear(self):
        return self.queue[-1]

    # Return the value of and remove the element passed to the method
    def pop(self, value):
        return self.queue.remove(value)

    # Return length of queue
    def length(self) -> int:
        return len(self.queue)

    # Return a copy of the queue
    def copy(self):
        return self.queue

# Class includes all specification for camera and methods used for FOV calculations
class Camera:
    def __init__(self, x, y, z, fov, focus, range):
        self.x = x
        self.y = y
        self.z = z
        self.fov = fov
        self.focus = focus
        self.min_range = range[0]
        self.max_range = range[1]
        self.checked_points = []

    # Return list of points not seen by the camera
    def check_fov(self, point, possible_obstructions, planar_equations, step, tC, tP):
        x, y, z = self.calculate_coordinates(point, step)
        confirmed_obstructions = self.check_for_obstructions(planar_equations, possible_obstructions, x, y, z, point, step, tC, tP)

        return confirmed_obstructions

    # Return array of coordinates for the ray between camera and point
    def calculate_coordinates(self, point, step):
        x = np.linspace(self.x, point[0], step, endpoint=False)
        y = np.linspace(self.y, point[1], step, endpoint=False)
        z = np.linspace(self.z, point[2], step, endpoint=False)

        return x, y, z

    # Checks if ray passes through any surface
    @staticmethod
    def check_for_obstructions(equations, surfaces, x, y, z, point,
                               step, threshold_crossing, threshold_point):
        for j in range(step):
            for k, l in enumerate(equations):
                planar_result = l[0] * (x[j] - l[1]) + l[2] * (y[j] - l[3]) + l[4] * (z[j] - l[5])
                crossing_plane = abs(planar_result) < threshold_crossing
                if crossing_plane:
                    # Calculates coordinates of margin-box around surface
                    corners_x = (np.min([surfaces[k][0][0], surfaces[k][1][0],
                                         surfaces[k][2][0], surfaces[k][3][0]]),
                                 np.max([surfaces[k][0][0], surfaces[k][1][0],
                                         surfaces[k][2][0], surfaces[k][3][0]]))
                    corners_y = (np.min([surfaces[k][0][1], surfaces[k][1][1],
                                         surfaces[k][2][1], surfaces[k][3][1]]),
                                 np.max([surfaces[k][0][1], surfaces[k][1][1],
                                         surfaces[k][2][1], surfaces[k][3][1]]))
                    corners_z = (np.min([surfaces[k][0][2], surfaces[k][1][2],
                                         surfaces[k][2][2], surfaces[k][3][2]]),
                                 np.max([surfaces[k][0][2], surfaces[k][1][2],
                                         surfaces[k][2][2], surfaces[k][3][2]]))
                    # Checks whether or not ray passes through surface and that the surface is not part of the same one
                    # as the point itself
                    confirmed_obstruction = (corners_x[0] - threshold_point) <= x[j] <= (corners_x[1] + threshold_point) and \
                             (corners_y[0] - threshold_point) <= y[j] <= (corners_y[1] + threshold_point) and \
                             (corners_z[0] - threshold_point) <= z[j] <= (corners_z[1] + threshold_point) and \
                             not (point[0] - threshold_point <= x[j] <= point[0] + threshold_point and
                                  point[1] - threshold_point <= y[j] <= point[1] + threshold_point and
                                  point[2] - threshold_point <= z[j] <= point[2] + threshold_point)

                    if confirmed_obstruction:
                        return True

        return False

    # Calculates field of view vectors and determines which points are within the field of view
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
                # Points that are not outside of field of view are within
                possible_points.append(i)

        return possible_points, points_outside_fov, unit_vectors

    # Return vector between camera and focus point
    def get_focus_vector(self):
        vector = np.array([self.focus[0]-self.x, self.focus[1]-self.y, self.focus[2]-self.z])
        unit_vector = vector / np.linalg.norm(vector)

        return unit_vector

    # Return vector perpendicular to the focus vector, along the XY-plane
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


# Class used to represent pallets as objects
class Pallet:
    def __init__(self, faces):
        self.faces = faces
        self.equations = get_planar_equations(faces)


# Shows points as seen or not seen in 3D-plot
def show_plots(seen_points, obstructed_points, cameras, unit_vectors=None, color=None):
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
    # Add points to plot
    if color is None:
        ax.scatter(X_seen, Y_seen, Z_seen, color="green", marker=".")
    else:
        ax.scatter(X_seen, Y_seen, Z_seen, color=color, marker=".")
    ax.scatter(X_obs, Y_obs, Z_obs, color="red", marker=".")
    # Add field of view vectors to the plot
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


# Calculate planar equations for each of the surfaces of the pallet
def get_planar_equations(surfaces):
    equations = []
    for j in surfaces:
        vector_1 = [j[1][0] - j[0][0], j[1][1] - j[0][1], j[1][2] - j[0][2]]
        vector_2 = [j[2][0] - j[0][0], j[2][1] - j[0][1], j[2][2] - j[0][2]]
        normal = np.cross(vector_1, vector_2)
        normal_hat = normal / (normal ** 2).sum() ** 0.5
        equations.append(get_plane(j[0], normal_hat))

    return equations


# Return the vectors representing the field of view of the camera
def get_camera_fov(camera):
    focus_vector = camera.get_focus_vector()
    focus_axis_1 = camera.get_focus_axis_1(focus_vector)
    focus_axis_2 = get_focus_axis_2(focus_vector, focus_axis_1)
    fov_vectors = get_fov_vectors(focus_axis_1, focus_axis_2, focus_vector, camera.fov)

    return fov_vectors

# Return components of planar equation as individual elements
def get_plane(point, normal):
    return normal[0], point[0], normal[1], point[1], normal[2], point[2]


# Return second focus axis perpendicular to the first and the focus vector
def get_focus_axis_2(focus_vector, focus_axis):
    return np.cross(focus_vector, focus_axis)


# Return field of view vectors for each camera as a list
def get_fov_vectors(vector_1, vector_2, focus_vector, fov):
    left_turn_radians = fov[0] * np.pi / 180 / 2
    top_turn_radians = fov[1] * np.pi / 180 / 2

    leftwards = np.tan(left_turn_radians)
    rightwards = -leftwards
    upwards = np.tan(top_turn_radians)
    downwards = -upwards

    vector_top_left = leftwards * vector_1 + upwards * vector_2 + focus_vector
    vector_bottom_left = leftwards * vector_1 + downwards * vector_2 + focus_vector
    vector_bottom_right = rightwards * vector_1 + downwards * vector_2 + focus_vector
    vector_top_right = rightwards * vector_1 + upwards * vector_2 + focus_vector

    return [vector_top_left, vector_bottom_left, vector_bottom_right, vector_top_right]



# Return list of seen points based on camera specifications and pallet surfaces
def check_fov(camera, current_camera, pallet, remaining_points, steps, threshold_point, threshold_surface, seen_points=[], queue=[], multi_processing=False):
    # If the setup specifies not using multi processing
    if not multi_processing:
        camera_progress = 0
        for j in remaining_points:
            point_obstructed = camera.check_fov(j, pallet.faces, pallet.equations, steps, threshold_surface, threshold_point)
            if point_obstructed:
                pass
            else:
                seen_points.append(j)
            camera_progress += 1
            print("Camera " + str(current_camera) + " " + str(
                round(camera_progress / len(remaining_points) * 100, 2)) + "% done")

        return seen_points
    # If setup specifies using multi processing
    else:
        while len(queue) > 0:
            if queue[0] not in camera.checked_points:
                current_point = queue.pop(0)
                camera.checked_points.append(current_point)
            elif len(queue) > 1:
                queue_copy = queue
                pos = 1
                rounds = 0
                checked = True
                while checked:
                    if pos < len(queue_copy):
                        if queue_copy[pos] not in camera.checked_points:
                            current_point = queue_copy[pos]
                            try:
                                queue.remove(current_point)
                            except ValueError:
                                pass
                            else:
                                camera.checked_points.append(current_point)
                                checked = False
                        else:
                            pos += 1
                    # When three passes of the remaining points do not return any unchecked points
                    elif rounds == 3:
                        print(f"Camera {current_camera} finished")
                        return

                    else:
                        pos = 1
                        rounds += 1
                        print(f"Camera {current_camera} finished {rounds} rounds without new points")
                        time.sleep(2)
            # If no points remain unseen
            else:
                print(f"Camera {current_camera} finished")
                return
            # Print remaining number of unseen points
            print(len(queue))
            point_obstructed = camera.check_fov(current_point, pallet.faces, pallet.equations, steps, threshold_surface, threshold_point)

            if point_obstructed:
                queue.append(current_point)
            else:
                seen_points.append(current_point)


# Import .txt-file and return surfaces and points as arrays
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
            # Every four points represent a surface
            for j in range(4):
                initial_points.append([float(pallet_surface_text[i * 12 + j * 3]),
                                        float(pallet_surface_text[i * 12 + 1 + j * 3]),
                                        float(pallet_surface_text[i * 12 + 2 + j * 3])])
                pallet_surface_temp.append([float(pallet_surface_text[i * 12 + j * 3]),
                                          float(pallet_surface_text[i * 12 + 1 + j * 3]),
                                          float(pallet_surface_text[i * 12 + 2 + j * 3])])
            pallet_surface_num.append(pallet_surface_temp)
        # Removes duplications of points
        for j in initial_points:
            if initial_points.count(j) > 1:
                initial_points.remove(j)
        # Add specified amount of points along sides of surfaces
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


# Setup and main program
def main():
    debug = False  # Debug mode
    fov_check = False  # FOV check mode
    use_multi_processing = False  # Multi processing mode
    steps = 5000  # Number of discretization points between camera and point
    additional_points = 10  # Additional points along each surface
    threshold_surface = 1/steps*20  # Margins of the planar equations
    threshold_point = 0.02  # Margin for deciding whether point is part of a crossed surface
    start_time = time.time()
    pallet_surface_num, initial_points = read_file("euro.txt", additional_points)
    cameras = []
    # Mode used to check if camera FOV is correct
    if debug:
        cameras.append(Camera(-10, -10, 0, (1, 1), (0, 0, 100), (0, 150)))
        p1 = Pallet(pallet_surface_num)

        points_in_fov, not_possible_points, unit_vectors = cameras[0].check_points_in_fov(initial_points)
        seen_points_camera = check_fov(cameras[0], 1, p1, points_in_fov, steps, threshold_point, threshold_surface)

        obstructed_points = []
        for i in initial_points:
            if i not in seen_points_camera:
                obstructed_points.append(i)
        show_plots(points_in_fov, not_possible_points, cameras, unit_vectors)
    # Mode used to check which points are within field of view. Does not check for obstructions
    elif fov_check:
        cameras.append(Camera(-80, -60, 0, (33, 25), (40, 60, 10), (50, 200)))
        points_in_fov, not_possible_points, unit_vectors = cameras[0].check_points_in_fov(initial_points)
        show_plots(points_in_fov, not_possible_points, cameras, unit_vectors, color="blue")
    # Checks all point to see whether or not they are seen by any of the cameras
    else:
        cameras.append(Camera(-38, -18, -44, (70, 55), (40, 40, 5), (25, 900)))
        cameras.append(Camera(-38, -18, 100, (70, 55), (40, 40, 5), (25, 900)))
        cameras.append(Camera(118, 138, -44, (70, 55), (40, 40, 5), (25, 900)))

        p1 = Pallet(pallet_surface_num)
        # Checks seen points using multi processing
        if use_multi_processing:
            manager = multiprocessing.Manager()
            seen_points = manager.list()
            queue = manager.list()

            for i in initial_points:
                queue.append(i)
            process_pool = []
            for i, c in enumerate(cameras):
                remaining_points, obstructed_points, unit_vectors = c.check_points_in_fov(initial_points)
                process_pool.append(multiprocessing.Process(target=check_fov, args=(c, i+1, p1, remaining_points, steps, threshold_point, threshold_surface, seen_points, queue, True)))
                process_pool[i].start()
            for process in process_pool:
                process.join()
        # Checks seen points one camera at a time
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
        # All unseen points is added to seperate list
        for i in initial_points:
            if i not in seen_points:
                obstructed_points.append(i)

        print("Total time used:", time.time() - start_time)
        # Displays results in 3D-plot
        try:
            show_plots(seen_points, obstructed_points, cameras)
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
