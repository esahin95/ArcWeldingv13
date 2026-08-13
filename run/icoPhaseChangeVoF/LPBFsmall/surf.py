import numpy as np
import matplotlib.pyplot as plt

import surf2stl as stl

M = np.loadtxt('postProcessing/traceSurface/0/traceSurface.dat', skiprows=2)


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.plot_trisurf(M[:, 0], M[:, 1], M[:, 2], linewidth=0.2, antialiased=True)
ax.set_aspect('equal', 'box')
ax.set_axis_off()

fig.tight_layout()

# export surface to a stl format file
X = M[:, 0].reshape(30,20)
Y = M[:, 1].reshape(30,20)
Z = M[:, 2].reshape(30,20)
stl.write('surf.stl', X, Y, Z)

plt.show()