#!/opt/local/bin/python2.7
##############################
# author Yuan Liu wendao@uw.edu

import os
import os.path
import sys
from math import pi,sin,cos,sqrt,fabs,exp,log
import numpy as np
from numpy import arccos,arctan2
from numpy import vstack,array
from numpy.random import random
from scipy.cluster.vq import kmeans,vq
from scipy.special import i0
import scipy.spatial.distance as sp
from pylab import *

# load the ref center
def load_centroid(fn, ncluster):
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
	return vstack(centroids)

# from 0 to 35 !! (0==36)
def pbc(bin):
	if bin>=36:
		bin=0
	if bin<0:
		bin=35
	return int(bin)

# [-180, 180) -> [0, 360) -> [0, 36)
# i.e.
# 0: (-185,-175)
# ...
# 35: (345,355)
def get_bb_bin_index(bb):
	bin1 = int((bb[0]+180)/10)
	bin2 = int((bb[1]+180)/10)
	return pbc(bin1), pbc(bin2)

#for debug
def show_bb_table(bbtable):
	for i in xrange(36):
		for j in xrange(36):
			print int(bbtable[i,j]), " ",
			#print bbtable[i,j], " ",
		print "\n",

#von Mises PDF
# dx in degree
def ksr(dx, kappa):
	dx = dx/180*pi
	return np.exp(kappa*np.cos(dx))/i0(kappa)/(2*pi)

#if the mean value of that angle close to the bundary
#dat point may shift by 2*pi which may cause the average wrong
def correct_pbc(dat):
	#check if the mean value around -180/180
	cosdat = np.cos(dat)
	meancos = np.average(cosdat)
	if meancos<-0.7:
		dat[dat>0] = dat[dat>0]-2*pi
	#or
	#dat[dat-acos(meancos)>pi] -2pi
	#dat[dat-acos(meancos)<-pi] +2pi

############################################################
##   main
############################################################
if len(sys.argv) < 3:
	print 'command line: script + AA + nrot'
	print 'kmeans.py <AA> <Nrot>'
	print 'pick up <Nrot> rot for <AA>'
	sys.exit()

#########
draw=False
#########

cutoff = 1.6
restype = sys.argv[1]
sigma = float(sys.argv[2])
scale = 2.0*sigma*sigma
ndata = int(sys.argv[3])

inpfile = "./split/" + restype + ".dat"

inp = open(inpfile, 'r')
lines = inp.readlines()
inp.close()

#cart coordinates for each sample
xyzlist = []
#internal coordinates for each sample
intlist = []
#psi, phi list
bblist = []

Nd = len(lines)
print Nd
prob = float(ndata)/Nd
print prob

for line in lines:
	pr = random()
	if pr>prob:
		continue
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

#prepare data (for clustering)
data = vstack(xyzlist)
save_data = data

Nd = len(data)
print Nd

print "affinity matrix"
Y = sp.squareform(sp.pdist(data, 'sqeuclidean'))
A = np.exp(-Y/scale)

print "D^(-1/2)"
D = np.zeros([Nd, Nd])
D_inv_sqr = np.zeros([Nd, Nd])
for i in xrange(Nd):
  D[i, i] = 0
  for j in xrange(Nd):
    D[i, i] += A[i, j]
  D_inv_sqr[i , i] = 1.0/sqrt(D[i, i])

L = ( D_inv_sqr.dot(A) ).dot( D_inv_sqr )
eig_vals, eig_vecs = linalg.eig(L)
eig_vecs = eig_vecs[eig_vals.argsort()]
eig_vals.sort()

x = range(20)
plot(x, eig_vals.real[-20:], 'ro')
plot(x, eig_vals.real[-20:], 'k--')
show()
