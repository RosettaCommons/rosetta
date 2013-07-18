import sys, os
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

inputfiles = sys.argv[1:]
xcol = 0
ycol = 2

fig = plt.figure()
i=1
for inputfile in inputfiles:
	inp = open(inputfile, 'r')
	lines = inp.readlines()
	inp.close()

	xdat = []
	ydat = []
	for line in lines:
		dat=line.split()
		x=float(dat[xcol])#*180/3.1416
		y=float(dat[ycol])*180/3.1416
		if x<9990 and y<9990:
			xdat.append(x)
			ydat.append(y)

#range = [[0, 180], [-180, 180]] #3-4
	range = [[0, 6], [-180, 180]] #2-4
#range = [[0, 6], [0, 180]] #2-3
	H, xedges, yedges = np.histogram2d(xdat, ydat, range=range, bins=[90, 180])
	extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]

	ctdatx = Counter(xdat)
	ctdaty = Counter(ydat)
	peakx = ctdatx.most_common(1)
	peaky = ctdaty.most_common(1)
	print inputfile, peakx[0][0], peaky[0][0]

	ax = fig.add_subplot(5,4,i)
	ax.imshow(H, extent=extent, interpolation='nearest')
#circ = plt.Circle((peaky[0][0], peakx[0][0]), radius=0.4, color='black')
#ax.add_patch(circ)
	
	ax.axes.get_xaxis().set_ticks([])
	ax.axes.get_yaxis().set_ticks([])
	ax.axes.set_xlabel(inputfile, fontsize=12)
#ax.set_aspect(15)#2-3
	ax.set_aspect(30)#2-4
	i = i+1

fig.savefig('2d.png')
plt.show()

