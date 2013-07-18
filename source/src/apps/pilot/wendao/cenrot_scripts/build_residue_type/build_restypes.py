#!/usr/bin/python
##############################
# author Yuan Liu wendao@uw.edu

import os
import sys

if len(sys.argv) == 1:
	print 'command line: script + AA'
	print 'build_restypes.py <AA>'
	print 'generate residue type file for <AA>'
	sys.exit()

restype = sys.argv[1]
print "#generating residue type file for ", restype

inpfile = "../centroid/residue_types/" + restype + ".params"
outfile = "residue_types/" + restype + ".params"

print "#template file " + inpfile
print "#output file " + outfile

#[0]:bond length
#[1]:bond angle (180-theta)
#[2]:diheral angle
ideal_dat_map = {
#"ALA":[0.367, 179.936, -0.015 ],
#"ARG":[3.903, 142.12, -117.351 ],
#"ASN":[1.872, 115.439, -66.191 ],
#"ASP":[1.849, 113.665, -70.402 ],
#"CYS":[1.695, 114.458, -68.313 ],
#"GLN":[2.436, 145.119, -65.48 ],
#"GLU":[2.416, 147.147, -69.455 ],
#"GLY":[0.654, 123.058, -179.998 ],
#"HIS":[2.633, 113.322, -58.006 ],
#"ILE":[1.304, 142.646, -88.003 ],
#"LEU":[1.878, 123.521, -79.62 ],
#"LYS":[3.159, 140.85, -69.85 ],
#"MET":[2.367, 150.073, -0.001 ],
#"PHE":[2.856, 114.371, -62.75 ],
#"PRO":[1.734, 87.496, -12.738 ],
#"SER":[1.249, 113.434, 63.152 ],
#"THR":[0.847, 129.494, -116.355 ],
#"TRP":[3.304, 112.333, -92.546 ],
#"TYR":[3.349, 113.671, -59.653 ],
#"VAL":[0.863, 128.776, -123.147 ],
"GLY":[0.000000,0.000000,0.000000],
"ALA":[0.000000,0.000000,0.000000],
"ARG":[3.980613,38.040012,-62.847846],
"ASN":[1.946896,67.287364,-67.236078],
"ASP":[1.934362,66.262497,-69.157053],
"CYS":[1.824346,67.931162,-78.859825],
"GLN":[2.521188,33.746074,-65.263501],
"GLU":[2.541004,33.574722,-65.144147],
"HIS":[2.668089,66.771385,-64.869587],
"ILE":[1.307011,39.071722,-84.333050],
"LEU":[1.923511,56.773401,-78.457443],
"LYS":[3.210392,39.915949,-67.452047],
"MET":[2.494268,31.951898,-71.173534],
"PHE":[2.883498,65.654597,-64.615458],
"PRO":[1.855479,115.601404,-7.007657],
"SER":[1.403274,67.935456,66.722850],
"THR":[0.850975,53.037438,12.512754],
"TRP":[3.296966,69.640641,-88.332341],
"TYR":[3.395916,66.019656,-63.744875],
"VAL":[0.869629,51.614101,-120.686781],
}

nbr_radius_map = {
#	"GLY": 2.636,
#	"ALA": 3.849,
#	"ARG": 3.974,
#	"ASN": 3.852,
#	"ASP": 3.867,
#	"CYS": 3.839,
#	"GLN": 3.842,
#	"GLU": 3.851,
#	"HIS": 3.842,
#	"ILE": 3.818,
#	"LEU": 3.841,
#	"LYS": 3.844,
#	"MET": 3.812,
#	"PHE": 3.836,
#	"PRO": 3.795,
#	"SER": 3.869,
#	"THR": 3.863,
#	"TRP": 3.834,
#	"TYR": 3.844,
#	"VAL": 3.860,
"GLY": 2.636,
"ALA": 3.794,
"ARG": 4.082,
"ASN": 3.782,
"ASP": 3.811,
"CYS": 3.778,
"GLN": 3.777,
"GLU": 3.798,
"HIS": 3.785,
"ILE": 3.765,
"LEU": 3.793,
"LYS": 3.794,
"MET": 3.772,
"PHE": 3.775,
"PRO": 3.794,
"SER": 3.793,
"THR": 3.782,
"TRP": 3.771,
"TYR": 3.780,
"VAL": 3.776,
}

inp = open(inpfile,'r')
lines = inp.readlines()
inp.close()

for line in lines:
	elem = line.split()
#check BOND CA CEN line
	if len(elem)==3 and elem[0]=="BOND" and elem[1]=="CA" and elem[2]=="CEN":
		if restype=="GLY":
			print "BOND CA CEN"
		else:
			print "BOND CB CEN"
		if restype!="GLY" and restype!="ALA":
			print "CHI 1 N CA CB CEN"
		continue
#check ICOOR_INT CEN line
	if len(elem)>2 and elem[0]=="ICOOR_INTERNAL" and elem[1]=="CEN":
		if restype=="GLY":
			print "ICOOR_INTERNAL CEN", ideal_dat_map[restype][2], ideal_dat_map[restype][1],ideal_dat_map[restype][0], "CA N C"
		else:
			print "ICOOR_INTERNAL CEN", ideal_dat_map[restype][2], ideal_dat_map[restype][1],ideal_dat_map[restype][0], "CB CA N"
		continue
#check NBR line
	if len(elem)>1 and elem[0]=="NBR_ATOM":
		if restype=="GLY":
			print "NBR_ATOM CA"
		else:
			print "NBR_ATOM CB"
		continue
	if len(elem)>1 and elem[0]=="NBR_RADIUS":
		print "NBR_RADIUS", nbr_radius_map[restype]
		continue
	print line,

