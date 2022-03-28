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
    ax.quiver(0, 0, 0, N_hat[0], N_hat[1], N_hat[2], color='r', linewidths=2)
    ax.quiver(0, 0, 0, AB_hat[0], AB_hat[1], AB_hat[2], color='b', linewidths=2)
    ax.quiver(0, 0, 0, AD_hat[0], AD_hat[1], AD_hat[2], color='b', linewidths=2)
    ax.plot_surface(x, y, z, color="purple")
    plt.show()


if __name__ == "__main__":
    main()