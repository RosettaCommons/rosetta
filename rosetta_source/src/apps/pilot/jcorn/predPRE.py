#!/bin/python
# predPRE.py - predict PRE and Iox/Ired from a PDB and labeling position
# Jacob Corn, 2/28/10
from math import *
from scipy import optimize
import sys

if( len(sys.argv) != 3 ):
	print "Usage: predPRE.py <pdbfilename> <labeled resnum + chain>. Eg - predPRE.py mypdb.pdb 125A"
	print "Output is in the format: resnum chain linewidth Iox Ired PRE"
	exit(1)

def ratio_func(x, R2, ratio):
	t=9e-3
	return R2*exp(-x*t) / (R2+x) - ratio

def fit( linewidth, ratio ):
	t=9e-3
	R2=linewidth/pi
	x0=40
	PRE = optimize.fsolve( ratio_func, x0, args=(R2, ratio))
	return PRE

K = 1.23e-32
larmor = 500e3
linewidth = 40 # temp hack
t=9e-3
R2 = linewidth/pi
nres=0

pdbname = sys.argv[1]
label = sys.argv[2]
label_chain = str(label[-1])
if( label_chain.isalpha() ):
	label_resnum = int(label[:-1])
else:
	label_chain = "A"
	label_resnum = int(label)


FILE=open( pdbname, 'r')
labeledx=0
labeledy=0
labeledz=0
for line in FILE.readlines():
	if( line=="\n"):
		continue
	fields = line.split()
	if( fields[0] != "ATOM" ):
		continue
	if( float(fields[5])==label_resnum and fields[4]==label_chain and fields[2]=="CA" ):
		labeledx=float(fields[6])
		labeledy=float(fields[7])
		labeledz=float(fields[8])
		nres += 1
	elif ( float(fields[5]) != label_resnum and fields[2]=="CA" ):
		nres += 1
if( labeledx==0 and labeledy==0 and labeledz==0 ):
	print "Residue %s not found!" % label
	exit(1)

FILE.seek(0) # return to the beginning of the file
MW = nres*110
kDa = MW/1000
tauc = 0.6e-9 * kDa

print label_resnum
for line in FILE.readlines():
	fields=line.split()
	if( fields[0] != "ATOM" ):
		continue
	if( fields[2]=="N" ):
		resnum=int(fields[5])
		chain=fields[4]
		if( resnum == label_resnum ):
			continue
		x=float(fields[6])
		y=float(fields[7])
		z=float(fields[8])
		dist=sqrt((x-labeledx)**2+(y-labeledy)**2+(z-labeledz)**2)
		dist = dist*10**-8

		Rconst = 4*tauc+((3*tauc)/(1+larmor**2 * tauc**2))
		numer = R2*exp(-t*(K/(dist**6)*Rconst))
		denom = R2+(K/(dist**6)*Rconst)
		ratio = float(numer)/denom
		PRE = fit( linewidth, ratio )
		Iox = ratio*100
		Ired = 100
		print "%s %.1f %.2f %.2f %.2f %s" % (resnum, linewidth, Iox, Ired, PRE, chain)

FILE.close()

