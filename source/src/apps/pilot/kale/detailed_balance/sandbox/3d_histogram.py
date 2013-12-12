#!/usr/bin/env python

import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('job')
parser.add_argument('residue')
parser.add_argument('--bins', '-b', type=int, default=360)
arguments = parser.parse_args()

from os import fork
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

fig = plt.figure()
data = np.load(os.path.join(arguments.job, 'histogram/torsions.npz'))

phi = data['phi' + arguments.residue]
psi = data['psi' + arguments.residue]

hist, xedges, yedges = np.histogram2d(phi, psi, bins=arguments.bins)

#xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:])
#xpos = xpos.flatten()/2.
#ypos = ypos.flatten()/2.
#zpos = np.zeros_like (xpos)
#dx = xedges [1] - xedges [0]
#dy = yedges [1] - yedges [0]
#dz = hist.flatten()
#ax = fig.add_subplot(111, projection='3d')
#ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')

plt.imshow(hist, interpolation='nearest') #, norm=mpl.colors.LogNorm())
plt.colorbar()
plt.xlabel('$\phi_%s$' % arguments.residue)
plt.ylabel('$\psi_%s$' % arguments.residue, rotation='horizontal')

ticks = np.arange(0, 360, 60)

if not fork(): plt.show()
