#!/usr/bin/env python

# This scripts demonstrates that randomly sampling a single torsion on each 
# move does not yield a uniform distribution of angles.  Sampling all of the 
# torsion angles does, however.  Because the non-uniformity is small and 
# inconsistent (i.e. there always is a non-uniformity, but it occurs in a 
# different places on different runs), I believe it is related to the fact that 
# each chosen angle is reported many times before it is changed again.  Maybe 
# this would amplify random deviations from the mean.
#
# Wait!  It totally would!  Because there are basically fewer samples than it 
# seems like there are.  If I run for algorithm for 15K iterations, only 15K 
# torsions are sampled.  But the KS test thinks there are 15K * N samples, 
# where N is the number of torsions in the peptide being simulated.  This 
# explains why the KS test keeps rejecting the null hypothesis.

import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('--iterations', '-i', type=int, default=15000)
parser.add_argument('--algorithm', '-a', default='one', choices=('one', 'all'))
arguments = parser.parse_args()

phi = [-56.9781, -57.051, -56.9634, -57.0267, -57.0011]
psi = [-47.0159, -46.9559, -46.9855, -46.9924, -47.0303]

header = 'apps.pilot.kale.monte_carlo:'
row = 5 * '%7.3f   '

def move_one_torsion():
    index = random.randrange(10)
    torsion = 360 * random.random()

    if index < 5:
        phi[index] = torsion;
    else:
        psi[index - 5] = torsion;

def move_all_torsions():
    for i in range(5):
        phi[i] = 360 * random.random()
        psi[i] = 360 * random.random()


for x in range(arguments.iterations):
    if arguments.algorithm == 'one':
        move_one_torsion()
    else:
        move_all_torsions()

    if x % 100 == 0:
        print header, 'phi     ', row % tuple(phi)
        print header, 'psi     ', row % tuple(psi)
        print header, 'random  ', 360 * random.random()


