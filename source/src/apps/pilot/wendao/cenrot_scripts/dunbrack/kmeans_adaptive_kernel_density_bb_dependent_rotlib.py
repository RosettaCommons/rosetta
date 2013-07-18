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
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
from scipy.special import i0

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm
import matplotlib.pyplot as plt

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
	bin1 = int((bb[0]+180)/10)
	bin2 = int((bb[1]+180)/10)
	return pbc(bin1), pbc(bin2)

def show_bb_table(bbtable):
	for i in xrange(36):
		for j in xrange(36):
			print int(bbtable[i,j]), " ",
			#print bbtable[i,j], " ",
		print "\n",

def ksr(dx, kappa):
	dx = dx/180*pi
	return np.exp(kappa*np.cos(dx))/i0(kappa)/(2*pi)

def correct_pbc(dat):
	#check if the mean value around -180/180
	cosdat = np.cos(dat)
	meancos = np.average(cosdat)
	if meancos<-0.7:
		dat[dat>0] = dat[dat>0]-2*pi
	#or
	#dat[dat-acos(meancos)>pi] -2pi
	#dat[dat-acos(meancos)<-pi] +2pi

##############################
# main
##############################
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
ncluster = int(sys.argv[2])
inpfile = "./split/" + restype + ".dat"

inp = open(inpfile, 'r')
lines = inp.readlines()
inp.close()

xyzlist = []
intlist = []
bblist = []

for line in lines:
	dats = line.split()
	dis = float(dats[0])
	ang = float(dats[1])
	dih = float(dats[2])
	psi = float(dats[3])
	phi = float(dats[4])
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

#assign
idx, dist = vq(save_data, centroids)

if draw:
	fig = plt.figure()

## Adaptive Kernel Desity Estimate
kappa_scale = 2.9
kappa = sqrt(savend)/kappa_scale #40~100 CYS~LEU
# this is a very arbitrary guess, may need to be change if change dataset

prhos = []
rhos = []
ncs = []
phis = []
psis = []
xphi = np.linspace(-180, 180, 36, endpoint=False)
ypsi = np.linspace(-180, 180, 36, endpoint=False)

################################
# !!! note: it's too slow, use bined rho as an approximation
#Step 0
#lambda_i for each data
# phis = np.zeros([savend])
# psis = np.zeros([savend])
# for nd in xrange(savend):
# 	#all data
# 	phis[nd] = bblist[nd][1]
# 	psis[nd] = bblist[nd][0]

# #for the residue type as a whole
# f_pp = np.zeros([savend])
# for nd in xrange(savend):
# 	phi = bblist[nd][1]
# 	psi = bblist[nd][0]
# 	dphis = phis - phi
# 	dpsis = psis - psi
# 	kphis = ksr(dphis,kappa)
# 	kpsis = ksr(dpsis,kappa)
# 	f_pp[nd] = np.dot(kphis, kpsis)/savend

# log_g = np.sum(np.log(f_pp))/savend
# g = exp(log_g)
# lambda_i = g / f_pp

#Step 1
#for each cluster, cal pilot density, N per degree^2 (rad?)
for i in xrange(ncluster):
	flag = idx==i #mark the data in cluster i
	n = sum(flag) #number of cluster i data
	ncs.append(n)
	philst = np.zeros([n])
	psilst = np.zeros([n])
	tn = 0
	for nd in xrange(savend):
		if(flag[nd]):
			philst[tn]=(bblist[nd][1])
			psilst[tn]=(bblist[nd][0])
			tn = tn+1
	assert (tn==n), "data number of cluster doesn't match"
	phis.append(philst)
	psis.append(psilst)
	rho = np.zeros([36,36])
	#calculate the contibution from all data points to each bin
	for x in xrange(36):
		for y in xrange(36):
			#cal for bin (x,y)
			phi = xphi[x]
			psi = ypsi[y]
			dphis = philst - phi
			dpsis = psilst - psi
			kphis = ksr(dphis,sqrt(n)/kappa_scale)#kappa=sqrt(n)/3
			kpsis = ksr(dpsis,sqrt(n)/kappa_scale)#kappa=sqrt(n)/3
			rho[x,y] = np.dot(kphis, kpsis)/n
	#pilot density
	prhos.append(rho)

#Step 2
rho_dat = np.zeros([savend])
all_phi = np.zeros([savend])
all_psi = np.zeros([savend])
#cal adaptive parameter for each data
for nd in xrange(savend):
	all_phi[nd] = bblist[nd][1]
	all_psi[nd] = bblist[nd][0]
	psi_bin, phi_bin = get_bb_bin_index(bblist[nd])
	for i in xrange(ncluster):
		rho_dat[nd]+=prhos[i][phi_bin,psi_bin]
log_rho = np.log(rho_dat)
log_g = np.sum(log_rho)/savend
#alpha=1/2 as paper suggest
lambda_i = np.exp( (log_g - log_rho)/2 )

#cross validation, opt kappa
#kappa = cross_validation_check(all_phi, all_psi, rho_dat)
#well, i am lazy, just choose by eye is fine, i guess

#lambda is the adaptive parameter (degree of smoothness)
#kappa for each data point in the whole dataset
kappa_i = kappa / lambda_i
#print lambda_i

#Step 3
#cal smoothed rho
rhos = []
for i in xrange(ncluster):
	rho = np.zeros([36,36])
	rhos.append(rho)

for x in xrange(36):
	for y in xrange(36):
		#cal for bin (x,y)
		phi = xphi[x]
		psi = ypsi[y]
		dphis = all_phi - phi
		dpsis = all_psi - psi
		kphis = ksr(dphis,kappa_i)
		kpsis = ksr(dpsis,kappa_i)
		#rho[x,y] = np.dot(kphis, kpsis)/n
		for i in xrange(ncluster):
			flag = idx==i
			rhos[i][x,y] = np.dot(kphis[flag], kpsis[flag])/ncs[i]

#draw
if draw:
	for i in xrange(ncluster):
		#draw it
		ax = fig.add_subplot(3,ncluster,i+1, projection='3d')
		x, y = np.meshgrid(xphi, ypsi)
		#print x.shape, y.shape, H.shape
		ax.set_aspect(1)
		ax.set_zlim3d(0, 5)   
		ax.view_init(50,45)
		#ax.zaxis.set_major_locator(LinearLocator(10))
		#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		surf=ax.plot_surface(x, y, rhos[i], rstride=1, cstride=1, 
			cmap=cm.coolwarm, linewidth=0, antialiased=False)
		#fig.colorbar(surf, shrink=0.5, aspect=5)

#Step 4:
#Bayes: P(r|phi, psi)
Pr = np.zeros([ncluster])
sumP = np.zeros([36,36])
for i in xrange(ncluster):
	Pr[i] = float(ncs[i])/sum(ncs)
	sumP = sumP + rhos[i]*Pr[i]

Pr_pp = []

for i in xrange(ncluster):
	#prpp = np.zeros([36,36])
	prpp = rhos[i]*Pr[i]/sumP
	Pr_pp.append(prpp)
	if draw:
		#draw it
		ax = fig.add_subplot(3,ncluster,i+ncluster+1, projection='3d')
		x, y = np.meshgrid(xphi, ypsi)
		#print x.shape, y.shape, H.shape
		ax.set_aspect(1)
		ax.view_init(50,45)
		ax.set_zlim3d(0, 1)
		#ax.zaxis.set_major_locator(LinearLocator(10))
		#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		surf=ax.plot_surface(x, y, prpp, rstride=1, cstride=1, 
			cmap=cm.coolwarm, linewidth=0, antialiased=False)

#Step 5:
#cal adaptive parameter for each bin
#and cal Mean({x}|phi,psi,rot), x=dis, ang, dih. and var
dis_lst=[]
ang_lst=[]
dih_lst=[]
vdis_lst=[]
vang_lst=[]
vdih_lst=[]
for i in xrange(ncluster):
	n = ncs[i]
	lambda_mat = np.zeros([36,36])

	log_g = 0
	for nd in xrange(n):
		phi_bin, psi_bin = get_bb_bin_index([phis[i][nd],psis[i][nd]])
		log_g += log(prhos[i][phi_bin, psi_bin])
	log_g /= n
	g = exp(log_g)
	#print g
	#cal lambda_phipsi
	for x in xrange(36):
		for y in xrange(36):
			lambda_mat[x,y] = (g/prhos[i][x,y])**(0.5)
	#specify kappa for each data
	kappa_lst = np.zeros([n])
	for nd in xrange(n):
		phi_bin, psi_bin = get_bb_bin_index([phis[i][nd],psis[i][nd]])
		kappa_lst[nd] = sqrt(n)/kappa_scale / lambda_mat[phi_bin, psi_bin] #kappa=sqrt(n)/3

	dis = np.zeros([n])
	ang = np.zeros([n])
	dih = np.zeros([n])
	#load data
	flag = idx==i
	tn = 0
	for nd in xrange(savend):
		if flag[nd]:
			dis[tn], ang[tn], dih[tn] = intlist[nd]
			tn += 1
	assert (tn==n), "data number of cluster doesn't match"

	#correct the dih data in case -180
	correct_pbc(dih)

	#calculate mean value for each bin mean=ave(x)
	dis_mat = np.zeros([36,36])
	ang_mat = np.zeros([36,36])
	dih_mat = np.zeros([36,36])
	#calculate vari(sigma^2) for each bin var=ave((mean-x)^2)
	vdis_mat = np.zeros([36,36])
	vang_mat = np.zeros([36,36])
	vdih_mat = np.zeros([36,36])
	for x in xrange(36):
		for y in xrange(36):
			phi = xphi[x]
			psi = ypsi[y]
			dphis = phis[i] - phi
			dpsis = psis[i] - psi
			kphis = ksr(dphis,kappa_lst)
			kpsis = ksr(dpsis,kappa_lst)
			kernel = kpsis*kphis
			denominator = np.sum(kernel)
			dis_mat[x,y] = np.dot(kernel, dis)/denominator
			ang_mat[x,y] = np.dot(kernel, ang)/denominator
			dih_mat[x,y] = np.dot(kernel, dih)/denominator
			#cal distance^2 between each data to mean
			#len_seq = distance * distance
			#var_mat[x,y] = np.dot(kernel, len_seq)/denominator
			ddis = dis - dis_mat[x,y]
			vdis_mat[x,y] = np.dot(kernel, ddis*ddis) / denominator
			dang = ang - ang_mat[x,y]
			vang_mat[x,y] = np.dot(kernel, dang*dang) / denominator
			ddih = dih - dih_mat[x,y]
			vdih_mat[x,y] = np.dot(kernel, ddih*ddih) / denominator

	#save the mean matrix for each cluster(rotamer)
	dis_lst.append(dis_mat)
	ang_lst.append(ang_mat)
	dih_lst.append(dih_mat)
	vdis_lst.append(vdis_mat)
	vang_lst.append(vang_mat)
	vdih_lst.append(vdih_mat)

if draw:
#draw
	for i in xrange(ncluster):
		#draw it
		ax = fig.add_subplot(3,ncluster,i+ncluster*2+1, projection='3d')
		x, y = np.meshgrid(xphi, ypsi)
		#print x.shape, y.shape, H.shape
		ax.set_aspect(1) 
		ax.view_init(50,45)
		#ax.zaxis.set_major_locator(LinearLocator(10))
		#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		surf=ax.plot_surface(x, y, dih_lst[i], rstride=1, cstride=1, 
			cmap=cm.coolwarm, linewidth=0, antialiased=False)
		#fig.colorbar(surf, shrink=0.5, aspect=5)
	
	plt.show()

#Step 6
#output
#rr2dd = (180/pi)*(180/pi)
rr2dd = 1.0
for phi_bin in xrange(36):
	for psi_bin in xrange(36):
		phi = int(xphi[phi_bin])
		psi = int(ypsi[psi_bin])
		for nc in xrange(ncluster):
			#pbc of dih
			dihval = dih_lst[nc][phi_bin, psi_bin]*180/pi
			if dihval < -180:
				dihval += 360
			elif dihval > 180:
				dihval -= 360
			print restype, phi, psi, nc+1, \
			Pr_pp[nc][phi_bin, psi_bin], \
			dis_lst[nc][phi_bin, psi_bin], \
			ang_lst[nc][phi_bin, psi_bin]*180/pi, \
			dihval, vdis_lst[nc][phi_bin, psi_bin], \
			vang_lst[nc][phi_bin, psi_bin]*rr2dd, \
			vdih_lst[nc][phi_bin, psi_bin]*rr2dd

