#!/usr/bin/env python

"""dbscan_adaptive_kernel_density_bb_dependent_rotlib.py

Usage:
	dakdbbd.py <AA> [<eps>] [<nbcut>]
	dakdbbd.py (-h | --help)
	dakdbbd.py --version

Options:
	-h --help    Show this screen.
	AA           residue type
	eps          eps parameter in dbscan
	nbcut        neighbour cutoff in dbscan
"""

from docopt import docopt
import numpy as np
from numpy import vstack
from math import pi,sin,cos,sqrt,fabs,exp,log
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from sklearn.cluster import DBSCAN

def parse_opt(opts):
	restyp = opts['<AA>']
	eps = 0.2
	if opts['<eps>']: eps = opts['<eps>']
	nbcut = 10
	if opts['<nbcut>']: nbcut = opts['<nbcut>']
	return  restyp, float(eps), int(nbcut)

def load_data(restyp):
	inpfile = "./split/" + restyp + ".dat"
	inp = open(inpfile, 'r')
	lines = inp.readlines()
	inp.close()

	xyzlist = []
	intlist = []
	bblist = []

	nl = len(lines)
	p = 1
	if nl>20000:
		#shink data
		p=(float(20000)/nl)
	print p

	for line in lines:
		if np.random.rand()>p: continue
		dats = line.split()
		dis = float(dats[0])
		ang = float(dats[1])
		dih = float(dats[2])
		psi = float(dats[3])
		phi = float(dats[4])
		#filter valid data
		if dis>0.1 and dis<9 and fabs(ang)<=180 and fabs(dih)<=180:
			if fabs(psi)<=180 and fabs(phi)<=180:
				ang = pi-ang
				dih = dih
				x = dis*sin(ang)*cos(dih)
				y = dis*sin(ang)*sin(dih)
				z = dis*cos(ang)
				xyzlist.append([x, y, z])
				intlist.append([dis, pi-ang, dih])
				bblist.append([psi, phi])

	ndat = len(xyzlist)
	if len(intlist)==ndat and len(bblist)==ndat:
		print ndat, "samples loaded from datafile", inpfile
	else:
		print "Wrong datafile", inpfile
		sys.exit(1)

	return vstack(xyzlist), vstack(intlist), vstack(bblist)

def plot_cluster(X, core_samples, labels, R):
	unique_labels = set(labels)
	colors = cm.Spectral(np.linspace(0, 1, len(unique_labels)))

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_aspect("equal")
	#k: cluster number -1 means outliers
	#col: color
	for k, col in zip(unique_labels, colors):
	    if k == -1:
	        # Black used for noise.
	        col = 'k'
	        markersize = 6
	    class_members = [index[0] for index in np.argwhere(labels == k)]
	    cluster_core_samples = [index for index in core_samples if labels[index] == k]

	    elem_ndx = 0
	    for index in class_members:
	    	elem_ndx += 1
	    	if elem_ndx % 10 != 0: continue
	        x = X[index]
	        if index in core_samples and k != -1:
	            markersize = 14
	        else:
	            markersize = 6
	        #pl.plot(x[0], x[1], 'o', markerfacecolor=col, markeredgecolor='k', markersize=markersize)
	        ax.scatter3D(x[0],x[1],x[2], c=col, s=markersize, marker='o')

	#draw a sphere frame
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	x = R*np.cos(u)*np.sin(v)
	y = R*np.sin(u)*np.sin(v)
	z = R*np.cos(v)
	ax.plot_wireframe(x, y, z, color="r")
	plt.show()

# main function
def main():
	opts = docopt(__doc__, version='0.1 Author: Yuan Liu wendao@uw.edu')
	restyp, eps, nbcut = parse_opt(opts)
	#print "AA =", restyp, " eps =", eps, " nbcut =", nbcut

	xyzdat, intdat, bbdat = load_data(restyp)
	#save_xyz = xyzdat

	#cluster!!
	db = DBSCAN(eps=eps, min_samples=nbcut).fit(xyzdat)
	core_samples = db.core_sample_indices_
	labels = db.labels_

	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	print('Estimated number of clusters: %d' % n_clusters_)

	plot_cluster(xyzdat, core_samples, labels, np.max(intdat[:,0]))

#main function
if __name__ == "__main__":	
	main()
