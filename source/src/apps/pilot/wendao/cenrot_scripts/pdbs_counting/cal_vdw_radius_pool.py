import multiprocessing, math
import os, sys, time
import numpy as np
import copy

nproc = 2

all_count_table = np.zeros((25,25,50))

def distance( r1, r2 ):
	dx = r1[0]-r2[0]
	dy = r1[1]-r2[1]
	dz = r1[2]-r2[2]
	return math.sqrt(dx*dx+dy*dy+dz*dz)

def get_pair_counts(pdb):
	binwidth = 0.2
	nbins = 50

	aa_to_id = { "N"   : 0,
        "CA"  : 1,
        "CB"  : 2,
        "C"   : 3,
        "O"   : 4,
        "ALA" : 5,
        "CYS" : 6,
        "ASP" : 7,
        "GLU" : 8,
        "PHE" : 9,
        "GLY" : 10,
        "HIS" : 11,
        "ILE" : 12,
        "LYS" : 13,
        "LEU" : 14,
        "MET" : 15,
        "ASN" : 16,
        "PRO" : 17,
        "GLN" : 18,
        "ARG" : 19,
        "SER" : 20,
        "THR" : 21,
        "VAL" : 22,
        "TRP" : 23,
        "TYR" : 24 }

	count_table = np.zeros((25,25,50))

	dummy_xyz = [ 0, 0, 0 ]
	dummy_res = [dummy_xyz,dummy_xyz,dummy_xyz,dummy_xyz,dummy_xyz,dummy_xyz]

	inp = open(pdb,'r')
	lines = inp.readlines()
	inp.close()
	
	#for one pdb
	coords = []
	centypes = []
	oldres = 0
	resnum = 0
	for line in lines:
		if line[:4] != "ATOM":
			continue
		#parse the ATOM line
		resndx = int(line[22:26])
		if resndx != oldres:
			#new res
			oldres = resndx
			resnum = resnum + 1
			#assert resndx==resnum
			coords.append(copy.deepcopy(dummy_res))
			restype = line[17:20]
			centypes.append(restype)

		X = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
		atom = line[13:16].strip()
		if atom == "N":
			coords[resnum-1][0]=X
		elif atom == "CA":
			coords[resnum-1][1]=X
		elif atom == "CB":
			coords[resnum-1][2]=X
		elif atom == "C":
			coords[resnum-1][3]=X
		elif atom == "O":
			coords[resnum-1][4]=X
		elif atom == "CEN":
			coords[resnum-1][5]=X

	for i in xrange(0, resnum):
		for j in xrange(i+8, resnum):
			res_list_i = [ 0, 1, 2, 3, 4, aa_to_id[centypes[i]] ]
			res_list_j = [ 0, 1, 2, 3, 4, aa_to_id[centypes[j]] ]
			for ai in xrange(0, 6):
				for aj in xrange(0, 6):
					dist = distance(coords[i][ai], coords[j][aj])
					distbin = int(dist/binwidth)
					if dist>0 and distbin<nbins and distbin>=0:
						count_table[res_list_i[ai]][res_list_j[aj]][distbin]=count_table[res_list_i[ai]][res_list_j[aj]][distbin]+1
						count_table[res_list_j[aj]][res_list_i[ai]][distbin]=count_table[res_list_j[aj]][res_list_i[ai]][distbin]+1
	return count_table

#main
pool = multiprocessing.Pool(processes=nproc)
ans = pool.map(get_pair_counts, sys.argv[1:])
#print all_count_table

for an in ans:
	all_count_table = all_count_table + an

for i in xrange(0,25):
	for j in xrange(0,25):
		buf = "cencen_counts(%d,%d,:) = [" % (i+1, j+1)
		for n in xrange(0,50):
			buf = buf + " %d" % all_count_table[i][j][n]
		buf = buf + "];"
		print buf
