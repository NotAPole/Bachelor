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
        focusVec = self.get_focus_vector()
        focusAxis1 = self.get_focus_axis(focusVec)
        focusAxis2 = get_focus_plane(focusVec, focusAxis1)
        fovVectors = get_FOV_vectors(focusAxis1, focusAxis2, focusVec, self.FOV)
        unitVectors = []
        for i in range(len(fovVectors)-1):
            normalVector = np.cross(fovVectors[i], fovVectors[i+1])
            unitVector = normalVector / (normalVector**2).sum()**0.5
            unitVectors.append(unitVector)
        normalVector = np.cross(fovVectors[-1], fovVectors[0])
        unitVector = normalVector / (normalVector ** 2).sum() ** 0.5
        unitVectors.append(unitVector)

        pointsOutsideFOV = []
        for i in points:
            dist = np.sqrt((i[0] - self.x) ** 2 + (i[1] - self.y)**2 + (i[2] - self.z) ** 2)
            if dist <= self.range:
                for j in unitVectors:
                    planar_result = np.dot(np.array([i[0], i[1], i[2]]) - np.array([self.x, self.y, self.z]), j)
                    if planar_result > 0:
                        pointsOutsideFOV.append(i)
                        break
                    else:
                        pass
            else:
                pointsOutsideFOV.append(i)

        possiblePoints = []
        for i in points:
            if i not in pointsOutsideFOV:
                possiblePoints.append(i)

        return possiblePoints, pointsOutsideFOV, unitVectors

    def get_focus_vector(self):
        vector = np.array([self.focus[0]-self.x, self.focus[1]-self.y, self.focus[2]-self.z])
        unitVector = vector / np.linalg.norm(vector)

        return unitVector

    @staticmethod
    def get_focus_axis(focusVector):
        spreadRad = np.pi
        rotationalMatrix = [
            [np.cos(spreadRad / 2), -np.sin(spreadRad / 2), 0],
            [np.sin(spreadRad / 2), np.cos(spreadRad / 2), 0],
            [0, 0, 0]
        ]

        vector = np.matmul(rotationalMatrix, focusVector)

        if (vector**2).sum()**0.5 == 0 and focusVector[2] == 1:
            return np.array([-1, 0, 0])

        elif (vector**2).sum()**0.5 == 0 and focusVector[2] == -1:
            return np.array([1, 0, 0])

        else:
            return vector / (vector**2).sum()**0.5

    @staticmethod
    def get_spread_vectors_z(focusVec, rotationPoints, spread): # NOT USED
        fourDimensionFocusVec = np.array([[focusVec[0]],
                                          [focusVec[1]],
                                          [focusVec[2]],
                                          [1]])
        spreadRad = spread*np.pi/180

        transposeMatrix1 = np.array([[1, 0, 0, -rotationPoints[0][0]],
                                     [0, 1, 0, -rotationPoints[0][1]],
                                     [0, 0, 1, -rotationPoints[0][2]],
                                     [0, 0, 0, 1]])
        invTransposeMatrix1 = np.linalg.inv(transposeMatrix1)

        rotationalVector = np.array([rotationPoints[0][0] - rotationPoints[1][0], rotationPoints[0][1] - rotationPoints[1][1], rotationPoints[0][2] - rotationPoints[1][2]])
        unitRotationalVector = rotationalVector / np.linalg.norm(rotationalVector)
        a, b, c = unitRotationalVector[0], unitRotationalVector[1], unitRotationalVector[2]
        d = np.sqrt(b**2 + c**2)
        rotationalMatrixX = np.array([[1, 0, 0, 0],
                                      [0, 1, 0, 0],
                                      [0, 0, 1, 0],
                                      [0, 0, 0, 1]])
        invRotationalMatrixX = np.linalg.inv(rotationalMatrixX)
        rotationalMatrixY = np.array([[d, 0, -a, 0],
                                      [0, 1, 0, 0],
                                      [a, 0, d, 0],
                                      [0, 0, 0, 1]])
        invRotationalMatrixY = np.linalg.inv(rotationalMatrixY)
        rotationalMatrixZ = np.array([[np.cos(spreadRad/2), np.sin(spreadRad/2), 0, 0],
                                      [-np.sin(spreadRad/2), np.cos(spreadRad/2), 0, 0],
                                      [0, 0, 1, 0],
                                      [0, 0, 0, 1]])
        rotatedVector = invTransposeMatrix1@invRotationalMatrixX@invRotationalMatrixY@rotationalMatrixZ@rotationalMatrixY@rotationalMatrixX@transposeMatrix1@fourDimensionFocusVec
        rotatedVector = [rotatedVector[0, 0], rotatedVector[1, 0], rotatedVector[2, 0]]
        unitRotatedVector = rotatedVector / np.linalg.norm(rotatedVector)

        return unitRotatedVector


class Pallet:
    def __init__(self, faces):
        self.faces = faces
        self.equations = get_planar_equations(faces)


def show_plots(seenPoints, obstructedPoints, cameras, unitVectors=None):
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
    plt.ylim(-100, 100)
    ax.set_zlim(-100, 100)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.scatter(X_seen, Y_seen, Z_seen, color="green", marker=".")
    ax.scatter(X_obs, Y_obs, Z_obs, color="red", marker=".")
    for i in cameras:
        ax.scatter(i.x, i.y, i.z, color="black")
        fovVectors = get_camera_FOV(i)
        ax.plot([i.x, i.x + fovVectors[0][0] * 75], [i.y, i.y + fovVectors[0][1] * 75],
                [i.z, i.z + fovVectors[0][2] * 75], color="blue")
        ax.plot([i.x, i.x + fovVectors[1][0] * 75], [i.y, i.y + fovVectors[1][1] * 75],
                [i.z, i.z + fovVectors[1][2] * 75], color="blue")
        ax.plot([i.x, i.x + fovVectors[2][0] * 75], [i.y, i.y + fovVectors[2][1] * 75],
                [i.z, i.z + fovVectors[2][2] * 75], color="blue")
        ax.plot([i.x, i.x + fovVectors[3][0] * 75], [i.y, i.y + fovVectors[3][1] * 75],
                [i.z, i.z + fovVectors[3][2] * 75], color="blue")
        if unitVectors is not None:
            for j in unitVectors:
                ax.plot([i.x, i.x + j[0] * 10], [i.y, i.y + j[1] * 10],
                        [i.z, i.z + j[2] * 10], color="orange")

    ax.view_init(30, 220)
    plt.show()


def get_planar_equations(surfaces):
    equations = []
    for j in surfaces:
        vector1 = [j[1][0] - j[0][0], j[1][1] - j[0][1], j[1][2] - j[0][2]]
        vector2 = [j[2][0] - j[0][0], j[2][1] - j[0][1], j[2][2] - j[0][2]]
        normal = np.cross(vector1, vector2)
        normal_hat = normal / (normal ** 2).sum() ** 0.5
        equations.append(get_plane(j[0], normal_hat))

    return equations


def get_camera_FOV(camera):
    focusVec = camera.get_focus_vector()
    focusAxis1 = camera.get_focus_axis(focusVec)
    focusAxis2 = get_focus_plane(focusVec, focusAxis1)
    fovVectors = get_FOV_vectors(focusAxis1, focusAxis2, focusVec, camera.FOV)

    return fovVectors


def get_plane(point, normal):
    return normal[0], point[0], normal[1], point[1], normal[2], point[2]


def get_focus_plane(focusVector, focusAxis):
    unitNormal = np.cross(focusVector, focusAxis)

    return unitNormal


def get_FOV_vectors(vector1, vector2, focusVector, FOV):
    leftTurnRad = FOV[0] * np.pi / 180 / 2
    topTurnRad = FOV[1] * np.pi / 180 / 2

    leftwards = np.sin(leftTurnRad) + (1-np.cos(leftTurnRad))/np.cos(leftTurnRad)
    rightwards = -leftwards
    upwards = np.sin(topTurnRad) + (1-np.cos(topTurnRad))/np.cos(topTurnRad)
    downwards = -upwards

    vectorTopLeft = leftwards * vector1 + upwards * vector2 + focusVector
    vectorBottomLeft = leftwards * vector1 + downwards * vector2 + focusVector
    vectorBottomRight = rightwards * vector1 + downwards * vector2 + focusVector
    vectorTopRight = rightwards * vector1 + upwards * vector2 + focusVector

    return [vectorTopLeft, vectorBottomLeft, vectorBottomRight, vectorTopRight]


def get_FOV_surfaces(point, vectors):
    surfaces = []
    for i in range(len(vectors)-1):
        surfaces.append([point, point + vectors[i], point + vectors[i+1]])
    surfaces.append([point, point + vectors[-1], point + vectors[0]])
    #print(surfaces)

    return surfaces


def get_spherical_coordinates(cartesian): # NOT USED
    r = np.sqrt(cartesian[0]**2 + cartesian[1]**2 + cartesian[2]**2)
    theta = np.arctan(cartesian[1]/cartesian[0])
    phi = np.arccos(cartesian[2]/r)

    return [r, theta, phi]


def get_cartesian_coordinates(sphere): # NOT USED
    cartesianVectors = []
    for i in sphere:
        x = i[0] * np.cos(i[1]) * np.sin(i[2])
        y = i[0] * np.sin(i[1]) * np.sin(i[2])
        z = i[0] * np.cos(i[2])
        cartesianVectors.append([x, y, z])

    return cartesianVectors


def get_rotated_vectors(sphere, FOV): # NOT USED
    thetaRotation = FOV[0]*np.pi/180/2
    phiRotation = FOV[1]*np.pi/180/2
    r, theta, phi = sphere
    topLeftVector = [r, theta+thetaRotation, phi-phiRotation]
    bottomLeftVector = [r, theta+thetaRotation, phi+phiRotation]
    bottomRightVector = [r, theta-thetaRotation, phi+phiRotation]
    topRightVector = [r, theta - thetaRotation, phi - phiRotation]

    return [topLeftVector, bottomLeftVector, bottomRightVector, topRightVector]


def main():
    debug = True
    useMultiProcessing = False
    steps = 1
    additionalPoints = 10
    thresholdSurface = 1/steps*20
    thresholdPoint = 0.02
    startTime = time.time()
    palletSurfaceNum, initialPoints = read_file("euro.txt", additionalPoints)
    if debug:
        cameras = []
        cameras.append(Camera(-10, -10, 0, (50, 36), (0, 0, 0), 150))
        p1 = Pallet(palletSurfaceNum)

        pointsInFov, notPossiblePoints, unitVectors = cameras[0].check_points_in_FOV(initialPoints)
        seenPointsCamera = check_FOV(cameras[0], 1, p1, pointsInFov, steps, thresholdPoint, thresholdSurface)

        obstructedPoints = []
        for i in initialPoints:
            if i not in seenPointsCamera:
                obstructedPoints.append(i)
        show_plots(pointsInFov, notPossiblePoints, cameras, unitVectors)


    else:
        cameras = []
        cameras.append(Camera(-74.7, -125.3, -86, (50, 25), (0, 0, 10), 260))
        #cameras.append(Camera(116.4, 257.2, -86, (70, 55), (0, 0, 0), 260))
        #cameras.append(Camera(60, 260, 94.4, (70, 55), (40, 60, 10), 260))
        #cameras.append(Camera(0, -50, 0, (70, 55), (0, 0, 20), 1000))
        p1 = Pallet(palletSurfaceNum)
        if useMultiProcessing:
            manager = multiprocessing.Manager()
            seenPoints = manager.list()
            processPool = []
            for i, c in enumerate(cameras):
                remainingPoints, obstructedPoints = c.check_points_in_FOV(initialPoints)
                processPool.append(multiprocessing.Process(target=check_FOV, args=(c, i+1, p1, remainingPoints, steps, thresholdPoint, thresholdSurface, seenPoints)))
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
                pointsInFov, _, _ = c.check_points_in_FOV(remainingPoints)
                seenPointsCamera = check_FOV(c, i+1, p1, pointsInFov, steps, thresholdPoint, thresholdSurface)
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
