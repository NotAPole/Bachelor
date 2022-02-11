"""
README

@author: Pål-André Furnes
this project is developed to be used for proof-of-concept in my bachelor thesis.
a pallet is represented by an array of 3D points.
a camera is represented by a single point in space facing a specified direction with a specified field of view.
the resulting 3D-plot show green and red points. green means seen, red means unseen.
"""
import time

import numpy
import numpy as np
import matplotlib.pyplot as plt


class Camera:
    def __init__(self, x, y, z, FOV, focus, range):
        self.x = x
        self.y = y
        self.z = z
        self.FOV = FOV
        self.focus = focus
        self.range = range

    def check_possible_obstructions(self, point, pallet_faces):
        posObs = []
        pointDistance = ((point[0] - self.x) ** 2 + (point[1] - self.y) ** 2 + (point[2] - self.z) ** 2) ** (1 / 2)
        for i in pallet_faces:
            posObs.append(i)
            pass
            for j in i:
                minimumDistance = ((j[0] - self.x) ** 2 + (j[1] - self.y) ** 2 + (j[2] - self.z) ** 2) ** (1 / 2)
                if minimumDistance <= pointDistance:
                    posObs.append(i)
                    break

        return posObs

    def check_FOV(self, point, possibleObstructions, equations, step, tC, tP):
        planarEquations = equations
        x, y, z = self.calculate_coordinates(point, step)
        conObs = self.check_for_obstructions(planarEquations, possibleObstructions, x, y, z, point, step, tC, tP)

        return conObs

    def calculate_coordinates(self, point, step):
        x = np.linspace(self.x, point[0], step)
        y = np.linspace(self.y, point[1], step)
        z = np.linspace(self.z, point[2], step)

        return x, y, z

    @staticmethod
    def check_for_obstructions(equations, surfaces, x, y, z, point, step, tC, tP):
        thresholdCrossing = tC
        thresholdPoint = tP
        for j in range(step):
            for k, l in enumerate(equations):
                planar_result = l[0] * (x[j] - l[1]) + l[2] * (y[j] - l[3]) + l[4] * (z[j] - l[5])
                crossingPlane = abs(planar_result) < thresholdCrossing
                if crossingPlane:
                    corners_x = (np.min([surfaces[k][0][0], surfaces[k][1][0], surfaces[k][2][0], surfaces[k][3][0]]),
                                 np.max([surfaces[k][0][0], surfaces[k][1][0], surfaces[k][2][0], surfaces[k][3][0]]))
                    corners_y = (np.min([surfaces[k][0][1], surfaces[k][1][1], surfaces[k][2][1], surfaces[k][3][1]]),
                                 np.max([surfaces[k][0][1], surfaces[k][1][1], surfaces[k][2][1], surfaces[k][3][1]]))
                    corners_z = (np.min([surfaces[k][0][2], surfaces[k][1][2], surfaces[k][2][2], surfaces[k][3][2]]),
                                 np.max([surfaces[k][0][2], surfaces[k][1][2], surfaces[k][2][2], surfaces[k][3][2]]))
                    conObs = (corners_x[0] - thresholdPoint) <= x[j] <= (corners_x[1] + thresholdPoint) and \
                             (corners_y[0] - thresholdPoint) <= y[j] <= (corners_y[1] + thresholdPoint) and \
                             (corners_z[0] - thresholdPoint) <= z[j] <= (corners_z[1] + thresholdPoint) and \
                             not (point[0] - thresholdPoint <= x[j] <= point[0] + thresholdPoint and
                                  point[1] - thresholdPoint <= y[j] <= point[1] + thresholdPoint and
                                  point[2] - thresholdPoint <= z[j] <= point[2] + thresholdPoint)

                    if conObs:
                        return True

        return False

    def check_points_in_FOV(self, points):
        possiblePoints = []
        pointsOutsideFOV = []
        unitFocusVector = self.get_focus_vector()

        xyVectors = self.get_spread_vectors_xy(unitFocusVector, self.FOV[0])
        #zVectors = self.get_spread_vectors_z(unitFocusVector, self.FOV[1])
        for i in points:
            dist = np.sqrt((i[0] - self.x) ** 2 + (i[1] - self.y)**2 + (i[2] - self.z) ** 2)
            if dist <= self.range:
                pass

        return possiblePoints, pointsOutsideFOV

    def get_focus_vector(self):
        vector = numpy.array([self.focus[0]-self.x, self.focus[1]-self.y, self.focus[2]-self.z])

        return vector / ((vector ** 2).sum() ** 0.5)

    @staticmethod
    def get_spread_vectors_xy(focusVector, spread):
        spreadRad = spread*np.pi/180
        print(focusVector, spread)
        rotationalMatrix1 = [
            [np.cos(spreadRad / 2), -np.sin(spreadRad / 2), 0],
            [np.sin(spreadRad / 2), np.cos(spreadRad / 2), 0],
            [0, 0, 1]
        ]
        rotationalMatrix2 = [
            [np.cos(-spread / 2), -np.sin(-spread / 2), 0],
            [np.sin(-spread / 2), np.cos(-spread / 2), 0],
            [0, 0, 1]
        ]

        vectors = [np.matmul(rotationalMatrix1, focusVector), np.matmul(rotationalMatrix2, focusVector)]

        return vectors

    @staticmethod
    def get_spread_vectors_z(focusVector, spread):
        spreadRad = spread*np.pi/180
        print(focusVector)
        return vectors

class Pallet:
    def __init__(self, faces):
        self.faces = faces
        self.equations = self.get_planar_equations(faces)

    @staticmethod
    def get_planar_equations(surfaces):
        equations = []
        for j in surfaces:
            vector1 = [j[1][0] - j[0][0], j[1][1] - j[0][1], j[1][2] - j[0][2]]
            vector2 = [j[2][0] - j[0][0], j[2][1] - j[0][1], j[2][2] - j[0][2]]
            normal = np.cross(vector1, vector2)
            normal_hat = normal / (normal ** 2).sum() ** 0.5
            equations.append(get_plane(j[0], normal_hat))

        return equations


def show_plots(seenPoints, obstructedPoints, cameras):
    X_seen, Y_seen, Z_seen = [], [], []
    for j in seenPoints:
        X_seen.append(j[0])
        Y_seen.append(j[1])
        Z_seen.append(j[2])

    X_obs, Y_obs, Z_obs = [], [], []
    for j in obstructedPoints:
        X_obs.append(j[0])
        Y_obs.append(j[1])
        Z_obs.append(j[2])

    plt.figure()
    ax = plt.axes(projection='3d')
    plt.xlim(-10, 100)
    plt.ylim(0, 150)
    ax.set_zlim(0, 20)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.scatter(X_seen, Y_seen, Z_seen, color="green", marker=".")
    ax.scatter(X_obs, Y_obs, Z_obs, color="red", marker=".")
    for i in cameras:
        ax.scatter(i.x, i.y, i.z, color="black")
    ax.view_init(30, 220)
    plt.show()


def get_plane(point, normal):
    return normal[0], point[0], normal[1], point[1], normal[2], point[2]


def main():
    steps = 1000
    thresholdSurface = 0.005
    thresholdPoint = 0.02
    startTime = time.time()
    palletSurfaceNum, remainingPoints = read_file("euro.txt")

    cameras = []
    cameras.append(Camera(40, 60, 50, (10, 10), (0, 0, 0), 150))
    # cameras.append(Camera(40, 60, -50, (10, 10), (0, 0, 0)))
    cameras.append(Camera(-20, -20, 7, (10, 10), (0, 0, 0), 150))
    p1 = Pallet(palletSurfaceNum)

    seenPoints = []
    currentCamera = 1
    for i in cameras:
        remainingPoints, obstructedPoints = i.check_points_in_FOV(remainingPoints)
        cameraProgress = 0
        for j in remainingPoints:
            pointObstructed = i.check_FOV(j, p1.faces, p1.equations, steps, thresholdSurface, thresholdPoint)
            if pointObstructed:
                obstructedPoints.append(j)
            else:
                seenPoints.append(j)
            cameraProgress += 1
            print("Camera " + str(currentCamera) + " " + str(
                round(cameraProgress / len(remainingPoints) * 100, 2)) + "% done")
        remainingPoints = obstructedPoints
        currentCamera += 1

    print("Total time used:", time.time() - startTime)

    show_plots(seenPoints, obstructedPoints, cameras)


def read_file(file):
    with open(file) as palletFile:
        palletSurfaceText = palletFile.read()
        palletSurfaceText = palletSurfaceText.replace("[", "")
        palletSurfaceText = palletSurfaceText.replace("]", "")
        palletSurfaceText = palletSurfaceText.replace(" ", "")
        palletSurfaceText = palletSurfaceText.split(",")
        palletSurfaceNum = []
        remainingPoints = []
        for i in range(int(len(palletSurfaceText) / 12.0)):
            palletSurfaceTemp = []
            for j in range(4):
                remainingPoints.append([float(palletSurfaceText[i * 12 + j * 3]),
                                        float(palletSurfaceText[i * 12 + 1 + j * 3]),
                                        float(palletSurfaceText[i * 12 + 2 + j * 3])])
                palletSurfaceTemp.append([float(palletSurfaceText[i * 12 + j * 3]),
                                          float(palletSurfaceText[i * 12 + 1 + j * 3]),
                                          float(palletSurfaceText[i * 12 + 2 + j * 3])])
            palletSurfaceNum.append(palletSurfaceTemp)
        for j in remainingPoints:
            if remainingPoints.count(j) > 1:
                remainingPoints.remove(j)
    return palletSurfaceNum, remainingPoints


if __name__ == "__main__":
    main()
