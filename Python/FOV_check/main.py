"""
README

@author: Pål-André Furnes
this project is developed to be used for proof-of-concept in a bachelor thesis.
a pallet is represented by an array of 3D points.
a camera is represented by a single point in space facing a specified direction and field of view.
the resulting plot show blue and red points. blue means seen, red means unseen.
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
        conObs = self.check_for_obstructions(planarEquations, possibleObstructions, x, y, z)

        return conObs

    def calculate_coordinates(self, point, step):
        x, y, z = [], [], []
        for j in range(step + 1):
            x.append((point[0] - self.x) / step * j + self.x)
            y.append((point[1] - self.y) / step * j + self.y)
            z.append((point[2] - self.z) / step * j + self.z)
            # print(x, y, z)

        return x, y, z

    def check_for_obstructions(self, equations, surfaces, x, y, z):
        threshold = 0.001
        for j in range(len(x)):
            for k, l in enumerate(equations):
                planar_result = l[0]*(x[j]-l[1]) + l[2]*(y[j]-l[3]) + l[4]*(z[j]-l[5])
                #print(planar_result)
                crossingPlane = abs(planar_result) < threshold
                if crossingPlane:
                    print("d")
                    corners_x = (surfaces[k][0][0], surfaces[k][2][0])
                    corners_y = (surfaces[k][0][1], surfaces[k][2][1])
                    corners_z = (surfaces[k][0][2], surfaces[k][2][2])
                    print(x[j], np.min(corners_x), np.max(corners_x))
                    conObs = np.min(corners_x) <= x[j] <= np.max(corners_x) and np.min(corners_y) <= y[j] <= np.max(corners_y) and np.min(corners_z) <= z[j] <= np.max(corners_z)
                    print(conObs)
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

    y_tilt = 50
    x_tilt = 24
    for j in range(2):
        plt.figure()
        ax = plt.axes(projection='3d')
        plt.xlim(0, 6)
        plt.ylim(0, 6)
        ax.set_zlim(0, 6)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.scatter(X_seen, Y_seen, Z_seen, color="green")
        ax.scatter(X_obs, Y_obs, Z_obs, color="red")
        for k in palletFaces:
            for l in k:
                ax.scatter(l[0], l[1], l[2], color="blue")
        ax.scatter(cameraPoint[0], cameraPoint[1], cameraPoint[2], color="black")
        ax.view_init(y_tilt, x_tilt)
        plt.show()
        plt.close()
        y_tilt -= 45
        x_tilt += 20

def get_planar_equations(surfaces):
    equations = []
    for j in surfaces:
        # print(i)
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
    seenPoints = []
    obstructedPoints = []

    c1 = Camera(0, 1, 0, (10, 10))
    p1 = Pallet([
        [[1, 0, 0], [1, 2, 0], [1, 2, 2], [1, 0, 2]]
    ])
    points = [(2, 1, 1), (2,1,0), (2,1,3), (2,1,4), (2,1,4.1)]
    for i in points:
        possibleObstructions = c1.check_possible_obstructions(i, p1.faces)
        pointObstructed = c1.check_FOV(i, possibleObstructions, 1000)
        print(pointObstructed)
        if pointObstructed:
            obstructedPoints.append(i)
        else:
            seenPoints.append(i)

    show_plots(seenPoints, obstructedPoints, p1.faces, (c1.x, c1.y, c1.z))
    print("Total time used:", time.time()-startTime)



