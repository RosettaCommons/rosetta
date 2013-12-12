#!/usr/bin/env python

import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('--bins', '-b', type=int, default=30)
arguments = parser.parse_args()

from os import fork
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

x, y = [], []

with open(arguments.input) as file:
    for line in file:
        phi, psi = line.split()
        x.append(float(phi))
        y.append(float(psi))

x = np.array(x)
y = np.array(y)

hist, xedges, yedges = np.histogram2d(x, y, bins=arguments.bins)
xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:])

hist /= np.max(hist)

xpos = xpos.flatten()/2.
ypos = ypos.flatten()/2.
zpos = np.zeros_like (xpos)

dx = xedges [1] - xedges [0]
dy = yedges [1] - yedges [0]
dz = hist.flatten()

#ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')

plt.imshow(hist, interpolation='nearest')
plt.colorbar()
plt.xlabel ('Phi0')
plt.ylabel ('Psi0')
plt.xlim(0, 360)
plt.ylim(0, 360)

if not fork(): plt.show()
