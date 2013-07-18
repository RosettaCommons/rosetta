#!/opt/local/bin/python2.7
##############################
# author Yuan Liu wendao@uw.edu

import os
import os.path
import sys
from math import pi,sin,cos,sqrt,fabs
import numpy as np
from numpy import arccos,arctan2
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq

#load the ref center
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

def pbc(bin):
	if bin>=36:
		bin=0
	if bin<0:
		bin=35
	return int(bin)

def get_bb_bin_index(bb):
	psi_bin = int((bb[0]+180)/10)
	phi_bin = int((bb[1]+180)/10)
	return pbc(psi_bin), pbc(phi_bin)

def show_bb_table(bbtable):
	for i in xrange(36):
		for j in xrange(36):
			print int(bbtable[i,j]), " ",
			#print bbtable[i,j], " ",
		print "\n",

#old way
def get_prior(rot_bb_data, n):
	#decompse
	marginal_psi = np.zeros([36])
	marginal_phi = np.zeros([36])
	for ncoli in xrange(36):
		for ncolj in xrange(36):
			marginal_psi[ncoli]+=rot_bb_data[ncoli,ncolj]
			marginal_phi[ncolj]+=rot_bb_data[ncoli,ncolj]

	#prior count "extra" for each bin
	extra = 5
	for ncol in xrange(36):
		marginal_psi[ncol]+=36*extra
		marginal_psi[ncol]/=(n*4+36*36*extra)
		marginal_phi[ncol]+=36*extra
		marginal_phi[ncol]/=(n*4+36*36*extra)

	return marginal_psi, marginal_phi, n*4+36*36

#in the dunbrack_1997 paper they use SVD to get the prior
#something wrong here, got big and negative number :(
def get_prior_svd(rot_bb_data, n):
	A = np.zeros([36*36, 36*2])
	b = np.zeros([36*36])
	for i in xrange(36):
		for j in xrange(36):
			A[j+i*36, i   ] = 1
			A[j+i*36, j+36] = 1
			b[j+i*36] = rot_bb_data[i,j]
	U, s, V = np.linalg.svd(A)
	c = np.dot(U.T,b) #1296*1
	y = np.zeros([36*2])
	for i in xrange(36*2):
		if s[i]<0.0001:
			y[i] = 0;
		else:
			y[i] = c[i]/s[i]
	#print y.shape
	#print V.shape
	xSVD = np.dot(V.T, y)
	#print xSVD
	marginal_psi = np.zeros([36])
	marginal_phi = np.zeros([36])
	for i in xrange(36):
		marginal_psi[i]=xSVD[i]
		marginal_phi[i]=xSVD[i+36]
	return marginal_psi, marginal_phi, n*4

if len(sys.argv) < 3:
	print 'command line: script + AA + nrot'
	print 'kmeans.py <AA> <Nrot>'
	print 'pick up <Nrot> rot for <AA>'
	sys.exit()

cutoff = 1.6
restype = sys.argv[1]
ncluster = int(sys.argv[2])
inpfile = "./split/" + restype + ".dat"

inp = open(inpfile, 'r')
lines = inp.readlines()
inp.close()

xyzlist = []
bblist = []

for line in lines:
	dats = line.split()
	dis = float(dats[0])
	ang = float(dats[1])
	dih = float(dats[2])
	psi = float(dats[3])
	phi = float(dats[4])
	if dis>0.1 and dis<9 and fabs(ang)<=180 and fabs(dih)<=180 and fabs(psi)<=180 and fabs(phi)<=180:
		ang = pi-ang
		dih = dih
		x = dis*sin(ang)*cos(dih)
		y = dis*sin(ang)*sin(dih)
		z = dis*cos(ang)
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
	centroids,_ = kmeans(data, ncluster)

idx, dist = vq(data, centroids)
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

# assign
idx, dist = vq(save_data, centroids)

# Bayesian
K0 = 0.5
center_int_lst = []
rot_bb_data_lst = []
for i in xrange(ncluster):
	#output the centroid
	x = centroids[i][0]
	y = centroids[i][1]
	z = centroids[i][2]
	r = sqrt(x*x+y*y+z*z)
	t = arccos(z/r)
	p = arctan2(y,x)

	flag = idx==i  #mark the data in cluster i
	n = sum(flag)  #number of cluster i data
	center_int_lst.append([r, 180-t*180/pi, p*180/pi]) #save the internal coordinates for the cluster

	rot_bb_data = np.zeros([36,36])
	for nd in xrange(savend):
		if (flag[nd]):
			#data in cluster i
			psi_bin, phi_bin = get_bb_bin_index(bblist[nd])
			rot_bb_data[pbc(psi_bin)  , pbc(phi_bin)  ] += 1
			rot_bb_data[pbc(psi_bin+1), pbc(phi_bin)  ] += 1
			rot_bb_data[pbc(psi_bin)  , pbc(phi_bin+1)] += 1
			rot_bb_data[pbc(psi_bin+1), pbc(phi_bin+1)] += 1

	#show_bb_table(rot_bb_data)

	marginal_psi, marginal_phi, newN = get_prior(rot_bb_data,n)

	for ncoli in xrange(36):
		for ncolj in xrange(36):
			rot_bb_data[ncoli,ncolj]+=marginal_psi[ncoli]*marginal_phi[ncolj]*newN

	rot_bb_data_lst.append(rot_bb_data)
	#print n
	#show_bb_table(rot_bb_data)
	#print restype, i+1, float(n)/float(savend), r, 180-t*180/pi, p*180/pi, np.std(dist[flag])

for ncoli in xrange(36):
	for ncolj in xrange(36):
		psi = int(ncolj*10 - 180)
		phi = int(ncoli*10 - 180)
		nrot = 0
		for nc in xrange(ncluster):
			nrot+=rot_bb_data_lst[nc][ncolj, ncoli]
		for nc in xrange(ncluster):
			print restype, phi, psi, nc+1, rot_bb_data_lst[nc][ncolj, ncoli]/nrot, center_int_lst[nc][0], center_int_lst[nc][1], center_int_lst[nc][2]

