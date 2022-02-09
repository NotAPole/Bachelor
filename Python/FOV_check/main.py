"""
README

@author: Pål-André Furnes
this project is developed to be used for proof-of-concept in my bachelor thesis.
a pallet is represented by an array of 3D points.
a camera is represented by a single point in space facing a specified direction with a specified field of view.
the resulting 3D-plot show green and red points. green means seen, red means unseen.
"""
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


class Camera:
    def __init__(self, x, y, z, FOV):
        self.x = x
        self.y = y
        self.z = z
        self.FOV = FOV

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

    def check_FOV(self, point, possibleObstructions, step):
        planarEquations = get_planar_equations(possibleObstructions)
        x, y, z = self.calculate_coordinates(point, step)
        conObs = self.check_for_obstructions(planarEquations, possibleObstructions, x, y, z, point, step)

        return conObs

    def calculate_coordinates(self, point, step):
        x = np.linspace(self.x, point[0], step)
        y = np.linspace(self.y, point[1], step)
        z = np.linspace(self.z, point[2], step)
        for i in range(len(x)):
            x[i] = x[i].round(2)
            y[i] = y[i].round(2)
            z[i] = z[i].round(2)
        return x, y, z

    def check_for_obstructions(self, equations, surfaces, x, y, z, point, step):
        thresholdCrossing = 0.008
        thresholdPoint = 0.5
        for j in range(step):
            for k, l in enumerate(equations):
                planar_result = l[0]*(x[j]-l[1]) + l[2]*(y[j]-l[3]) + l[4]*(z[j]-l[5])
                crossingPlane = abs(planar_result) < thresholdCrossing
                if crossingPlane:
                    corners_x = (np.min([surfaces[k][0][0], surfaces[k][1][0], surfaces[k][2][0], surfaces[k][3][0]]), np.max([surfaces[k][0][0], surfaces[k][1][0], surfaces[k][2][0], surfaces[k][3][0]]))
                    corners_y = (np.min([surfaces[k][0][1], surfaces[k][1][1], surfaces[k][2][1], surfaces[k][3][1]]), np.max([surfaces[k][0][1], surfaces[k][1][1], surfaces[k][2][1], surfaces[k][3][1]]))
                    corners_z = (np.min([surfaces[k][0][2], surfaces[k][1][2], surfaces[k][2][2], surfaces[k][3][2]]), np.max([surfaces[k][0][2], surfaces[k][1][2], surfaces[k][2][2], surfaces[k][3][2]]))
                    #if not(point[0]-thresholdPoint <= x[j] <= point[0]+thresholdPoint and point[1] - thresholdPoint <= y[j] <= point[1] + thresholdPoint and point[2] - thresholdPoint <= z[j] <= point[2] + thresholdPoint):
                        #print([x[j], y[j], z[j]], point, surfaces[k][0], surfaces[k][1], surfaces[k][2], surfaces[k][3],
                        #  (corners_x[0] - 0.1) <= x[j] <= (corners_x[1] + thresholdPoint),
                        #  (corners_y[0] - thresholdPoint) <= y[j] <= (corners_y[1] + thresholdPoint),
                        #  (corners_z[0] - thresholdPoint) <= z[j] <= (corners_z[1] + thresholdPoint))
                    conObs = (corners_x[0] - 0.1) <= x[j] <= (corners_x[1] + thresholdPoint) and (corners_y[0] - thresholdPoint) <= y[j] <= (corners_y[1] + thresholdPoint) and (corners_z[0] - thresholdPoint) <= z[j] <= (corners_z[1] + thresholdPoint) and not(point[0]-thresholdPoint <= x[j] <= point[0]+thresholdPoint and point[1] - thresholdPoint <= y[j] <= point[1] + thresholdPoint and point[2] - thresholdPoint <= z[j] <= point[2] + thresholdPoint)

                    if conObs:
                        return True

        return False


class Pallet:
    def __init__(self, faces):
        self.faces = faces


def show_plots(seenPoints, obstructedPoints, palletFaces, cameraPoint):
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
    ax.scatter(cameraPoint[0], cameraPoint[1], cameraPoint[2], color="black")
    ax.view_init(30, 220)
    plt.show()


def get_planar_equations(surfaces):
    equations = []
    for j in surfaces:
        # print(j)
        vector1 = [j[1][0] - j[0][0], j[1][1] - j[0][1], j[1][2] - j[0][2]]
        vector2 = [j[2][0] - j[0][0], j[2][1] - j[0][1], j[2][2] - j[0][2]]
        # print(vector1, vector2)
        normal = np.cross(vector1, vector2)
        normal_hat = normal / (normal ** 2).sum() ** 0.5
        # print(normal_hat)
        equations.append(get_plane(j[0], normal_hat))

    return equations


def get_plane(point, normal):
    return normal[0], point[0], normal[1], point[1], normal[2], point[2]


if __name__ == "__main__":
    startTime = time.time()
    with open("euro.txt") as palletFile:
        palletSurfaceText = palletFile.read()
        palletSurfaceText = palletSurfaceText.replace("[", "")
        palletSurfaceText = palletSurfaceText.replace("]", "")
        palletSurfaceText = palletSurfaceText.replace(" ", "")
        palletSurfaceText = palletSurfaceText.split(",")
        palletSurfaceNum = []
        points = []
        for i in range(int(len(palletSurfaceText)/12.0)):
            palletSurfaceTemp = []
            for j in range(4):
                points.append([float(palletSurfaceText[i*12+j*3]), float(palletSurfaceText[i*12+1+j*3]), float(palletSurfaceText[i*12+2+j*3])])
                palletSurfaceTemp.append([float(palletSurfaceText[i*12+j*3]), float(palletSurfaceText[i*12+1+j*3]), float(palletSurfaceText[i*12+2+j*3])])
            palletSurfaceNum.append(palletSurfaceTemp)

    seenPoints = []
    obstructedPoints = []

    c1 = Camera(62.5, 60, 50, (10, 10))
    p1 = Pallet(palletSurfaceNum)
    #show_plots(points, obstructedPoints, p1.faces, (c1.x, c1.y, c1.z))
    #point = points[-5:-4]
    #for i in range(0, len(points), 4):
    #    print([points[i], points[i+1], points[i+2], points[i+3]])
    #    show_plots([points[i], points[i+1], points[i+2], points[i+3]], [], p1.faces, (c1.x, c1.y, c1.z))

    for i in points:
        if points.count(i) > 1:
            points.remove(i)

    for i in points:
        #possibleObstructions = c1.check_possible_obstructions(i, p1.faces)
        pointObstructed = c1.check_FOV(i, p1.faces, 3000)
        if pointObstructed:
            obstructedPoints.append(i)
        else:
            seenPoints.append(i)

    #show_plots(seenPoints, points[:-40], p1.faces, (c1.x, c1.y, c1.z))
    print("Total time used:", time.time() - startTime)
    show_plots(seenPoints, obstructedPoints, p1.faces, (c1.x, c1.y, c1.z))




