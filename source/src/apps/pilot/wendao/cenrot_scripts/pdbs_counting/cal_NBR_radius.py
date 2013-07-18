#!/usr/bin/python
##############################
# author Yuan Liu wendao@uw.edu

import os
import sys
import numpy as np
from math import pi
argv=sys.argv

if len(argv) == 1:
	print 'command line: script + PDBs'
	print 'cal_NBR_radius.py <struct.pdb>'
	print 'calculate the NBR_radius for each sidechain'
	sys.exit()

inputfiles = argv[1:]
masslst = {'C':12.0107, 'O':15.9994, 'N':14.00674, 'S':32.066, 'H':1.00794}

def parse_res_lines(lines, nres):
	restype = lines[0][17:20].strip()
	ar = np.array([0,0,0])   #CA
	br = np.array([0,0,0])   #CB
	cnr = np.array([0,0,0])  #CEN
	bbnr = np.array([0,0,0]) #N
	bbor = np.array([0,0,0]) #O
	bbcr = np.array([0,0,0]) #C
	hr = np.array([0,0,0])   #H
	nm=0; am=0; bm=0; hm=0; cnm=0; bom=0; bcm=0;

	for l in lines:
		elem = l[76:78].strip()
		atom = l[12:16].strip()
		if len(elem)==1:
			mass = masslst[elem]
		else:
			mass = 1.0
		x = float(l[30:38])
		y = float(l[38:46])
		z = float(l[46:54])

		if atom == "N":
			bbnr = np.array([x,y,z])
			nm = mass
		elif atom == "CA":
			ar = np.array([x,y,z])
			am = mass
		elif atom == "CB":
			br = np.array([x,y,z])
			bm = mass
		elif atom == "H":
			hr = np.array([x,y,z])
			hm = mass
		elif atom == "CEN":
			cnr = np.array([x,y,z])
			cnm = mass
		elif atom == "C":
			bbcr = np.array([x,y,z])
			bcm = mass
		elif atom == "O":
			bbor = np.array([x,y,z])
			bom = mass

	nbr_r=0
	if restype == "GLY":
		nbr = ar
		nbrm = am
	else:
		nbr = br
		nbrm = bm

	if nbrm>0:
		if nm >0:
			d = np.linalg.norm(bbnr-nbr)
			nbr_r = max(nbr_r, d)
		if hm >0:
			d = np.linalg.norm(hr-nbr)
			nbr_r = max(nbr_r, d)
		if am >0:
			d = np.linalg.norm(ar-nbr)
			nbr_r = max(nbr_r, d)
		if cnm >0:
			d = np.linalg.norm(nbr-cnr)
			nbr_r = max(nbr_r, d)
		if bom >0:
			d = np.linalg.norm(bbor-nbr)
			nbr_r = max(nbr_r, d)
		if bcm >0:
			d = np.linalg.norm(bbcr-nbr)
			nbr_r = max(nbr_r, d)
		
	print restype, nbr_r
	return

for inputfile in inputfiles:
	inp = open(inputfile,'r')
	lines = inp.readlines()
	inp.close()

	#the offset should be 0
	current_res = 1
	reslines = []
	for line in lines:
		if line[:4]=="ATOM":
			res_num = int(line[22:26])
			if res_num == current_res:
				reslines.append(line)
			else:
				#parse the old one
				parse_res_lines(reslines, current_res)
				reslines = []
				reslines.append(line)
				current_res = res_num
	#the last one
	parse_res_lines(reslines, res_num)

