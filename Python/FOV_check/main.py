"""
README

@author: Pål-André Furnes
This project is developed to be used for proof-of-concept in my bachelor thesis.
A pallet is represented by an array of 3D points.
A camera is represented by a single point in space facing a specified direction with a specified field of view.
The resulting 3D-plot show green and red points. green means seen, red means unseen.
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing


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
                #print(planar_result, k, l)
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
                        #print("obstructed", (x[j], y[j], z[j]), point)
                        return True

        return False

    def check_points_in_FOV(self, points):
        possiblePoints = []
        pointsOutsideFOV = []
        unitFocusVector = self.get_focus_vector()

        xyVectors = self.get_spread_vectors_xy(unitFocusVector, self.FOV[0])
        #plt.figure()
        #ax = plt.axes(projection='3d')
        #plt.xlim(-10, 100)
        #plt.ylim(0, 150)
        #ax.set_zlim(0, 20)
        #ax.plot((xyVectors[0][0]+self.x, xyVectors[0][0]*100+self.x), (xyVectors[0][1]+self.y, xyVectors[0][1]*100+self.y), (xyVectors[0][2]+self.z, xyVectors[0][2]*100+self.z))
        #plt.show()
        #zVectors = self.get_spread_vectors_z(unitFocusVector, self.FOV[1])
        for i in points:
            dist = np.sqrt((i[0] - self.x) ** 2 + (i[1] - self.y)**2 + (i[2] - self.z) ** 2)
            if dist <= self.range:
                xLimits = [np.min([xyVectors[0][0] * dist, xyVectors[1][0] * dist]), np.max([xyVectors[0][0] * dist, xyVectors[1][0] * dist])]
                yLimits = [np.min([xyVectors[0][1] * dist, xyVectors[1][1] * dist]), np.max([xyVectors[0][1] * dist, xyVectors[1][1] * dist])]

                if xLimits[0] <= i[0]+self.x <= xLimits[1] and yLimits[0] <= i[1]+self.y <= yLimits[1]:
                    print("inside")
                    possiblePoints.append(i)
                else:
                    print(unitFocusVector)
                    print(xyVectors)
                    print(yLimits, i[1])
                    pointsOutsideFOV.append(i)
            else:
                pointsOutsideFOV.append(i)

        return possiblePoints, pointsOutsideFOV

    def get_focus_vector(self):
        vector = np.array([self.focus[0]-self.x, self.focus[1]-self.y, self.focus[2]-self.z])

        return vector / ((vector ** 2).sum() ** 0.5)

    @staticmethod
    def get_spread_vectors_xy(focusVector, spread):
        spreadRad = spread*np.pi/180
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
        vectors = []
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
    plt.xlim(-100, 100)
    plt.ylim(-100, 200)
    ax.set_zlim(-100, 50)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.scatter(X_seen, Y_seen, Z_seen, color="green", marker=".")
    ax.scatter(X_obs, Y_obs, Z_obs, color="red", marker=".")
    for i in cameras:
        ax.scatter(i.x, i.y, i.z, color="black")
    ax.view_init(30, 220)
    plt.show()


def show_camera_FOV(camera):
    pass


def get_plane(point, normal):
    return normal[0], point[0], normal[1], point[1], normal[2], point[2]


def main():
    debug = False
    useMultiProcessing = False
    steps = 10000
    additionalPoints = 10
    thresholdSurface = 1/steps*20
    thresholdPoint = 0.02
    startTime = time.time()
    palletSurfaceNum, initialPoints = read_file("euro.txt", additionalPoints)

    if debug:
        camera = Camera(-45, -45, 0, (10, 10), (0, 0, 0), 9000)
        show_camera_FOV(camera)

    else:
        cameras = []
        cameras.append(Camera(-74.7, -125.3, -86, (10, 10), (0, 0, 0), 9000))
        cameras.append(Camera(116.4, 257.2, -86, (10, 10), (0, 0, 0), 9000))
        cameras.append(Camera(60, 260, 94.4, (10, 10), (0, 0, 0), 9000))
        #cameras.append(Camera(40, 60, -100, (10, 10), (0, 0, 0), 9000))
        p1 = Pallet(palletSurfaceNum)
        if useMultiProcessing:
            manager = multiprocessing.Manager()
            seenPoints = manager.list()
            processPool = []
            for i, c in enumerate(cameras):
                #remainingPoints, obstructedPoints = c.check_points_in_FOV(initialPoints)
                processPool.append(multiprocessing.Process(target=check_FOV, args=(c, i+1, p1, initialPoints, steps, thresholdPoint, thresholdSurface, seenPoints)))
                processPool[i].start()
            for process in processPool:
                process.join()

            for i in seenPoints:
                if seenPoints.count(i) > 1:
                    seenPoints.remove(i)
        else:
            seenPoints = []
            remainingPoints = initialPoints
            for i, c in enumerate(cameras):
                seenPointsCamera = check_FOV(c, i+1, p1, remainingPoints, steps, thresholdPoint, thresholdSurface)
                remainingPoints = []
                for j in seenPointsCamera:
                    seenPoints.append(j)
                for j in initialPoints:
                    if j not in seenPoints:
                        remainingPoints.append(j)

        obstructedPoints = []
        for i in initialPoints:
            if i not in seenPoints:
                obstructedPoints.append(i)

        print("Total time used:", time.time() - startTime)

        show_plots(seenPoints, obstructedPoints, cameras)


def check_FOV(c, currentCamera, pallet, remainingPoints, steps, thresholdPoint, thresholdSurface, seenPoints=[]):
    cameraProgress = 0
    for j in remainingPoints:
        pointObstructed = c.check_FOV(j, pallet.faces, pallet.equations, steps, thresholdSurface, thresholdPoint)
        if pointObstructed:
            pass
        else:
            seenPoints.append(j)
        cameraProgress += 1
        print("Camera " + str(currentCamera) + " " + str(
            round(cameraProgress / len(remainingPoints) * 100, 2)) + "% done")

    return seenPoints


def read_file(file, additionalPoints):
    with open(file) as palletFile:
        palletSurfaceText = palletFile.read()
        palletSurfaceText = palletSurfaceText.replace("[", "")
        palletSurfaceText = palletSurfaceText.replace("]", "")
        palletSurfaceText = palletSurfaceText.replace(" ", "")
        palletSurfaceText = palletSurfaceText.split(",")
        palletSurfaceNum = []
        initialPoints = []
        for i in range(int(len(palletSurfaceText) / 12.0)):
            palletSurfaceTemp = []
            for j in range(4):
                initialPoints.append([float(palletSurfaceText[i * 12 + j * 3]),
                                        float(palletSurfaceText[i * 12 + 1 + j * 3]),
                                        float(palletSurfaceText[i * 12 + 2 + j * 3])])
                palletSurfaceTemp.append([float(palletSurfaceText[i * 12 + j * 3]),
                                          float(palletSurfaceText[i * 12 + 1 + j * 3]),
                                          float(palletSurfaceText[i * 12 + 2 + j * 3])])
            palletSurfaceNum.append(palletSurfaceTemp)
        for j in initialPoints:
            if initialPoints.count(j) > 1:
                initialPoints.remove(j)

        for j in palletSurfaceNum:
            if j[0][0] != j[1][0]:
                for k in np.linspace(j[0][0], j[1][0], additionalPoints, endpoint=False):
                    initialPoints.append([k, j[0][1], j[0][2]])
                for k in np.linspace(j[2][0], j[3][0], additionalPoints, endpoint=False):
                    initialPoints.append([k, j[2][1], j[2][2]])

                if j[1][1] != j[2][1]:
                    for k in np.linspace(j[1][1], j[2][1], additionalPoints, endpoint=False):
                        initialPoints.append([j[1][0], k, j[1][2]])
                    for k in np.linspace(j[3][1], j[0][1], additionalPoints, endpoint=False):
                        initialPoints.append([j[3][0], k, j[3][2]])

                if j[1][2] != j[2][2]:
                    for k in np.linspace(j[1][2], j[2][2], additionalPoints, endpoint=False):
                        initialPoints.append([j[1][0], j[1][1], k])
                    for k in np.linspace(j[3][2], j[0][2], additionalPoints, endpoint=False):
                        initialPoints.append([j[3][0], j[3][1], k])

            if j[0][1] != j[1][1]:
                for k in np.linspace(j[0][1], j[1][1], additionalPoints, endpoint=False):
                    initialPoints.append([j[0][0], k, j[0][2]])
                for k in np.linspace(j[2][1], j[3][1], additionalPoints, endpoint=False):
                    initialPoints.append([j[3][0], k, j[2][2]])

                if j[1][0] != j[2][0]:
                    for k in np.linspace(j[1][0], j[2][0], additionalPoints, endpoint=False):
                        initialPoints.append([k, j[1][1], j[1][2]])
                    for k in np.linspace(j[3][0], j[0][0], additionalPoints, endpoint=False):
                        initialPoints.append([k, j[3][1], j[3][2]])

                if j[1][2] != j[2][2]:
                    for k in np.linspace(j[1][2], j[2][2], additionalPoints, endpoint=False):
                        initialPoints.append([j[1][0], j[1][1], k])
                    for k in np.linspace(j[3][2], j[0][2], additionalPoints, endpoint=False):
                        initialPoints.append([j[3][0], j[3][1], k])

            if j[0][2] != j[1][2]:
                for k in np.linspace(j[0][2], j[1][2], additionalPoints, endpoint=False):
                    initialPoints.append([j[0][0], j[0][1], k])
                for k in np.linspace(j[2][0], j[3][0], additionalPoints, endpoint=False):
                    initialPoints.append([[j[2][0]], j[2][1], k])

                if j[1][0] != j[2][0]:
                    for k in np.linspace(j[1][0], j[2][0], additionalPoints, endpoint=False):
                        initialPoints.append([k, j[1][1], j[1][2]])
                    for k in np.linspace(j[3][0], j[0][0], additionalPoints, endpoint=False):
                        initialPoints.append([k, j[3][1], j[3][2]])

                if j[1][1] != j[2][1]:
                    for k in np.linspace(j[1][1], j[2][1], additionalPoints, endpoint=False):
                        initialPoints.append([j[1][0], k, j[1][2]])
                    for k in np.linspace(j[3][1], j[0][1], additionalPoints, endpoint=False):
                        initialPoints.append([j[3][0], k, j[3][2]])

        for j in initialPoints:
            if initialPoints.count(j) > 1:
                initialPoints.remove(j)

    return palletSurfaceNum, initialPoints


if __name__ == "__main__":
    main()
