import matplotlib.pyplot as plt
import numpy as np


def main():
    x = np.outer(np.linspace(0, 2, 10), np.ones(10))
    y = x.copy().T  # transpose
    z = y
    A = np.array([0, 0, 0])
    B = np.array([4, 0, 0])
    C = np.array([4, 4, 4])
    D = np.array([0, 4, 4])
    AB = B-A
    AD = D-A
    N = np.cross(AB, AD)
    AB_hat = AB / (AB ** 2).sum() ** 0.5
    AD_hat = AD / (AD ** 2).sum() ** 0.5
    N_hat = N / (N ** 2).sum() ** 0.5
    plt.figure()
    ax = plt.axes(projection='3d')
    plt.xlim(-1.01, 1.01)
    plt.ylim(-1.01, 1.01)
    ax.set_zlim(-1.01, 1.01)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    focus = np.array([0, np.sqrt(0.5), np.sqrt(0.5)])
    vector_2 = np.cross(focus, np.array([-1, 0, 0]))
    #ax.quiver(0, 0, 0, N_hat[0], N_hat[1], N_hat[2], color='r', linewidths=2)
    ax.quiver(0, 0, 0, focus[0], focus[1], focus[2], color='r', linewidths=2)
    ax.quiver(focus[0], focus[1], focus[2], focus[0]-1, 0, 0, color='b', linewidths=2)
    ax.quiver(focus[0], focus[1], focus[2], vector_2[0], vector_2[1], vector_2[2], color='b', linewidths=2)
    #ax.plot_surface(x, y, z, color="purple")
    plt.show()


if __name__ == "__main__":
    main()