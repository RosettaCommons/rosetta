#!/usr/bin/python
# Jacob Corn, Feb 2010.
# PRE equations from Battiste and Wagner, 2000
# requires scipy (http://www.scipy.org/)
from math import *
from scipy import optimize
import sys

def ratio_func(x, R2, ratio):
	t=9e-3
	return R2*exp(-x*t) / (R2+x) - ratio

def fit( linewidth, ratio ):
	t=9e-3
	R2=linewidth/pi
	x0=40
	PRE = optimize.fsolve( ratio_func, x0, args=(R2, ratio))
	return PRE

if( len(sys.argv) != 3 ):
	print "Usage: PREcst.py <file> <MW in Da>"
	print "<file> format:"
	print "labeledres1"
	print "crosspeak1 linewidth Iox Ired"
	print "crosspeak2 linewidth Iox Ired"
	print "labeledres2"
	print "crosspeak1 linewidth Iox Ired"
	exit(1)

MW = float(sys.argv[2])
kDa = MW / 1000
tauc = 0.6e-9 * kDa
K = 1.23e-32
larmor = 500e3


fname = sys.argv[1]
FILE = open(fname, "r")
for line in FILE.readlines():
	if( line=="\n"):
		continue
	param = line.split(" ")
	if( len(param) == 1 ):
		labeledres = int(param[0])
		continue
	peak = int(param[0])
	linewidth = float(param[1])
	Iox = float(param[2])
	Ired = float(param[3])
	ratio = Iox/Ired

	# large ratios made into negative constraints (see Battiste and Wagner)
	# lb = 30, ub = 1000, sd=1
	if ratio > 0.85:
		#print peak, ratio, PRE, 1000
		ratio = 0.85
		PRE = fit( linewidth, ratio )
		Rconst = 4*tauc+((3*tauc)/(1+larmor**2 * tauc**2))
		r = (K/PRE*Rconst)**(float(1)/6)*10**8

		print "AtomPair CA %d N %d BOUNDED %0.1f 1000.0 1.0 0.5 tag" % (labeledres, peak, r)
		continue

	# very low ratios are counted as 0.01 and made into a bounded func
	# lb=0, ub=r(0.01), sd=1
	if ratio <= 0.01:
		ratio = 0.01
		PRE = fit( linewidth, ratio )
		Rconst = 4*tauc+((3*tauc)/(1+larmor**2 * tauc**2))
		r = (K/PRE*Rconst)**(float(1)/6)*10**8

		print "AtomPair CA %d N %d BOUNDED 0.0 %0.1f 1.0 0.5 tag" % (labeledres, peak, r)
		continue

	PRE = fit( linewidth, ratio )

	Rconst = 4*tauc+((3*tauc)/(1+larmor**2 * tauc**2))
	r = (K/PRE*Rconst)**(float(1)/6)*10**8
	#print peak, ratio, PRE, 1000
	# well-defined ratios are made into harmonic func with sd 4
	print "AtomPair CA %d N %d HARMONIC %0.1f 4.0" % (labeledres, peak, r)
FILE.close()

