#!/opt/local/bin/python2.7
##############################
# author Yuan Liu wendao@uw.edu

import os
import os.path
import sys
from math import pi,sin,cos,sqrt
import numpy as np
from numpy import arccos,arctan2
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from collections import Counter


def get_common_int_coords_2(x, y, z):
	nd = len(x)
	rdat = []
	tdat = []
	pdat = []
	for i in xrange(nd):
		rdat.append(sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]))
		tdat.append(arccos(z[i]/rdat[i]))
		pdat.append(arctan2(y[i],x[i]))
	crdat = Counter(rdat)
	ctdat = Counter(tdat)
	cpdat = Counter(pdat)
	peakr = crdat.most_common(1)
	peakt = ctdat.most_common(1)
	peakp = cpdat.most_common(1)
	return [peakr[0][0], peakt[0][0], peakp[0][0]]
#return [np.median(rdat), np.median(tdat), np.median(pdat)]

def get_common_int_coords(x, y, z):
	nd = len(x)
	cxdat = Counter(x)
	cydat = Counter(y)
	czdat = Counter(z)
	peakx = cxdat.most_common(1)
	peaky = cydat.most_common(1)
	peakz = czdat.most_common(1)
	x0 = peakx[0][0];
	y0 = peaky[0][0];
	z0 = peakz[0][0];
	r = sqrt(x0*x0+y0*y0+z0*z0)
	t = arccos(z0/r)
	p = arctan2(z0,x0)
	return [r, t, p]

def load_centroid(fn,ncluster):
	inp = open(fn,'r')
	lines = inp.readlines()
	inp.close()
	centroids = []
	for line in lines:
		dats = line.split()
		dis = float(dats[3])
		ang = (180-float(dats[4]))*pi/180
		dih = float(dats[5])*pi/180
		x = dis*sin(ang)*cos(dih)
		y = dis*sin(ang)*sin(dih)
		z = dis*cos(ang)
		centroids.append([x,y,z])
	while len(centroids) < ncluster:
		centroids.append([0,-1,-2.5])
	return vstack(centroids)

if len(sys.argv) < 3:
	print 'command line: script + AA + nrot'
	print 'kmeans.py <AA> <Nrot>'
	print 'pick up <Nrot> rot for <AA>'
	sys.exit()

cutoff = 1.6
restype = sys.argv[1]
ncluster = int(sys.argv[2])
inpfile = "./split/" + restype + ".dat"

inp = open(inpfile,'r')
lines = inp.readlines()
inp.close()
skip = int(len(lines)/2000)

colorlist = ['r', 'b', 'g', 'y', 'c', 'm', 'k', 'w', 'r']
marklist = ['o', '^', '+', '>', (5,2), (5,0)]
#intlist = []
xyzlist = []
bblist = []

for line in lines:
	dats = line.split()
	dis = float(dats[0])
	ang = float(dats[1])
	dih = float(dats[2])
	psi = float(dats[3])
	phi = float(dats[4])
	if dis>0.1 and dis<9990 and ang<9990 and dih<9990:
		ang = pi-ang
		dih = dih
		x = dis*sin(ang)*cos(dih)
		y = dis*sin(ang)*sin(dih)
		z = dis*cos(ang)
		#intlist.append([dis, ang, dih])
		xyzlist.append([x,y,z])
		bblist.append([psi, phi])

#pre data
data = vstack(xyzlist)
save_data = data
#idat = vstack(intlist)
#cluster
fn = "ref/" + restype
if os.path.isfile(fn):
	centroids = load_centroid(fn, ncluster)
	centroids,_ = kmeans(data, centroids)
else:
	centroids,_ = kmeans(data,ncluster)

idx, dist = vq(data,centroids)
nd = len(data)
savend = nd

## cluster, find the centroid
while True: 
	data = data[dist<cutoff]
	#idat = idat[dist<cutoff]
	centroids,_ = kmeans(data, centroids)
	idx, dist = vq(data,centroids)
	if nd <= len(data):
		nd = len(data)
		break
	else:
		nd = len(data)
	if cutoff>0.6: cutoff=cutoff-0.2

## assign
idx, dist = vq(save_data,centroids)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in xrange(ncluster):
	#plot
	c = colorlist[i]
	#m = marklist[i]
	m = '.'

	flag = idx==i

	xs = save_data[flag, 0]
	ys = save_data[flag, 1]
	zs = save_data[flag, 2]
	n = len(xs)
	ax.scatter(xs[0:n:skip], ys[0:n:skip], zs[0:n:skip], c=c)#,marker=m
	ax.scatter3D(centroids[i,0],centroids[i,1],centroids[i,2], c='black', s=20, marker='s')

	#output the centroid
	x = centroids[i][0]
	y = centroids[i][1]
	z = centroids[i][2]
	r = sqrt(x*x+y*y+z*z)
	t = arccos(z/r)
	p = arctan2(y,x)
	print restype, i+1, float(n)/float(savend), r, 180-t*180/pi, p*180/pi, np.std(dist[idx==i])

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()

