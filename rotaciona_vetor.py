# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 19:29:43 2022

@author: Marlipe
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


a = np.array([0, 1, 0])

theta = 90

rotation_radians = np.radians(theta)
rotation_axis = np.array([0, 0, 1])
rotation_vector = rotation_radians * rotation_axis
rotation = R.from_rotvec(rotation_vector)
rotated_vec = rotation.apply(a)

phi = 135

b = rotation_axis
c = rotated_vec

rotation_radians = np.radians(phi)
rotation_axis = np.array(a)
rotation_vector = rotation_radians * rotation_axis
rotation = R.from_rotvec(rotation_vector)
rotated_vec = rotation.apply(c)

d = rotated_vec

plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Acústica de Raios - Sala de Aula')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')

ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])

ax.plot([0, a[0]], [0, a[1]], [0, a[2]],'-b')
ax.plot([0, b[0]], [0, b[1]], [0, b[2]],'--r')
ax.plot([0, c[0]], [0, c[1]], [0, c[2]],'-.g')
ax.plot([0, d[0]], [0, d[1]], [0, d[2]],':m')






