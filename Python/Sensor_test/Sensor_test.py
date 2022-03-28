import serial
import time
from math import sin, cos, radians
import numpy as np
import matplotlib

matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ser = serial.Serial('COM3', 9600, timeout=1)
time.sleep(2)
run = True
counter = 0

mode = 'quat'

# definerer først fire hjørner og lager så fire overflater med disse hjørnene
corners = [[-1, 0.5, 0.5], [+1, 0.5, 0.5], [0, -0.5, 0.5], [0, 0.5, -0.5]]
face1 = [corners[0], corners[1], corners[2]]
face2 = [corners[0], corners[1], corners[3]]
face3 = [corners[0], corners[2], corners[3]]
face4 = [corners[1], corners[2], corners[3]]

# samler alle overflatene i ett numpy array
vertices = np.array(face1 + face2 + face3 + face4, dtype=float)

# lager 3d objektet
ob = Poly3DCollection(vertices, linewidths=1, alpha=0.2)

ob.set_facecolor([0.5, 0.5, 1])
ob.set_edgecolor([0, 0, 0])
plt.ion()
fig = plt.figure()
plt.pause(1)
# Legger til subplot til figuren
ax = fig.add_subplot(111, projection='3d')
plt.show()

# legger 3d objektet til subplottet
ax.add_collection3d(ob)

s = [-2, -2, -2, 2, 2, 2]
ax.auto_scale_xyz(s, s, s)
plt.pause(1)

message = True
while run:
    line = ser.readline()
    try:
        decoded_string = line.decode()
    except:
        print('Serial data not found')
        message = False
    if message:
        decoded_string = decoded_string.strip()
        if mode == 'euler':
            if 'Orientation:' in decoded_string:
                angles_str = decoded_string.replace('Orientation: ', '')
                angles_str = angles_str.split(',')
                angles = []
                for i in range(len(angles_str)):
                    angles.append(float(angles_str[i]))

                x = radians(angles[0])
                y = radians(angles[1])
                z = radians(angles[2])
                '''a1 = [-sin(z) * sin(x) + cos(y) * cos(x) * cos(z), sin(z) * cos(x) + cos(y) * sin(x) * cos(z),
                      -cos(z) * sin(y)]
                a2 = [-cos(z) * sin(x) - cos(y) * cos(x) * sin(z), cos(z) * cos(x) - cos(y) * sin(x) * sin(z),
                      sin(z) * sin(y)]
                a3 = [sin(y) * cos(x), sin(y) * cos(x), cos(y)]'''
                a1 = [cos(z) * cos(y), cos(z) * sin(y)*sin(x) - sin(z) * cos(x),
                      cos(z) * sin(y) * cos(x) + sin(z)+sin(x)]
                a2 = [sin(z) * cos(y), sin(z) * sin(y) * sin(x) + cos(z) * cos(x),
                      sin(z) * sin(y) * cos(x) - cos(z) * sin(x)]
                a3 = [-sin(y), cos(y) * sin(x), cos(y) * cos(x)]
                rot_mat = [a1, a2, a3]

                new_vertices = np.matmul(vertices, rot_mat)
                ob.set_verts(new_vertices)
                counter += 1
                plt.pause(0.01)

        if mode == 'quat':
            if 'Quaternion:' in decoded_string:
                quat_str = decoded_string.replace('Quaternion: ', '')
                quat_str = quat_str.split(',')
                quat = []
                for i in range(len(quat_str)):
                    quat.append(float(quat_str[i]))
                q0 = quat[0]
                q1 = quat[1]
                q2 = quat[2]
                q3 = quat[3]
                print(quat)
                a1 = [2 * (q0 ** 2 + q1 ** 2) - 1, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)]
                a2 = [2 * (q1 * q2 + q0 * q3), 2 * (q0 ** 2 + q2 ** 2) - 1, 2 * (q2 * q3 - q0 * q1)]
                a3 = [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), 2 * (q0 ** 2 + q3 ** 2) - 1]
                rot_mat = [a1, a2, a3]

                new_vertices = np.matmul(vertices, rot_mat)
                ob.set_verts(new_vertices)
                counter += 1
                plt.pause(0.01)
    message = True
    if counter == 1000:
        run = False
ser.close()