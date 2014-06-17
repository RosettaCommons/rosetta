#!/usr/bin/env python

import sys

fin_names=sys.argv[1:]

for fin_name in fin_names:
	fin = open(fin_name)
	if (fin_name[-4:] == ".pdb"):
		ext = "pdb"
	else:
		print "Unknown file extension"
		exit
	fout_name = fin_name[:-4] + "_val2thr" + fin_name[-4:]
	fout = open(fout_name, 'w')
	on_val = False
	for line in fin.readlines() :
		if line[:3] not in ["ATO","TER","END"]:
			continue
		res = line[17:20]
		atom = line[12:16]
		if res == "VAL":
			if (not on_val):
				on_val = True
				print "swapping val for thr"
			res = "THR"

			if atom in ["HG12","HG13"]:
				continue
			elif atom == "HG11":
				atom = "HG1"
			elif atom == " CG1":
				atom = " OG1"
		else:
			on_val = False
			
		newline = line[:12] + atom + " " + res + line[20:]
		fout.write(newline)
	fout.close()
	fin.close()
