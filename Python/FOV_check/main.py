"""
README

this project is developed by Pål-André Furnes to be used for proof-of-concept in a bachelor thesis.
a pallet is represented by an array of 3D points.
a camera is represented by a single point in space facing a specified direction and field of view.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


class camera:
    def __init__(self, x, y, z, FOV):
        self.x = x
        self.y = y
        self.z = z
        self.FOV = FOV

    def check_relevant_obstructions(self, point, pallet_faces):
        possibleObstructions = []
        pointDistance = ((point[0] - self.x) ** 2 + (point[1] - self.y) ** 2 + (point[2] - self.z) ** 2) ** (1 / 2)
        for i in pallet_faces:
            for j in i:
                minimumDistance = ((j[0] - self.x) ** 2 + (j[1] - self.y) ** 2 + (j[2] - self.z) ** 2) ** (1 / 2)
                if minimumDistance <= pointDistance:
                    possibleObstructions.append(i)
        print(possibleObstructions)


class pallet:
    def __init__(self, faces):
        self.faces = faces


if __name__ == "__main__":
    c1 = camera(1, 0, 0, (10, 10))
    p1 = pallet([[
        [3, 3, 3], [5, 5, 5], [6, 6, 6], [4, 4, 4]
    ]])
    # c1.check_relevant_obstructions((3, 3, 3), p1.faces)

    plt.figure()

    # Make data.
    X = [3, 3, 5, 5]
    Y = [3, 5, 5, 3]
    Z = [0, 0, 0, 0]

    # Plot the surface.
    ax = plt.axes(projection='3d')
    plt.xlim(0, 6)
    plt.ylim(0, 6)
    ax.set_zlim(0,6)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.scatter(X, Y, Z)
    ax.view_init(45, 180)
    plt.show()
    plt.show()
    plt.close()

    plt.figure()
    ax = plt.axes(projection='3d')
    plt.xlim(0, 6)
    plt.ylim(0, 6)
    ax.set_zlim(0,6)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.scatter(X, Y, Z)
    ax.view_init(45, 280)
    plt.show()
