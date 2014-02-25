#!/usr/bin/env python
# encoding: utf-8

import os, sys
import argparse
import biophysics
import numpy, matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Liberation Serif'
mpl.rcParams['axes.labelsize'] = 16

from pylab import *

parser = argparse.ArgumentParser()
parser.add_argument('job')
parser.add_argument('--limit', '-l', type=int)
parser.add_argument('--quiet', '-q', action='store_true')
parser.add_argument('--solutions', '-s', action='store_true')
parser.add_argument('--output', '-o')
parser.add_argument('--dpi', type=int, default=300)
arguments = parser.parse_args()

def backbone_filter(model, atom, position):
    first_residue = model.atoms[0]['residue-id']
    real_pivot = first_residue + third_pivot - 1

    if atom['residue-id'] == real_pivot:
        return atom['atom-name'] == 'CA'
    else:
        return False

def distance(coordinates):
    return numpy.linalg.norm(reference - coordinates)


# Calculate the how far the atoms beyond that pivot moved.
directory = os.path.join(arguments.job, 'trajectory')
paths = os.listdir(directory)
paths.sort()

if not paths:
    print "No output found in '%s'." % directory
    raise SystemExit
if arguments.limit:
    paths = paths[:arguments.limit + 1]

count = len(paths)
models = arange(count)
error = zeros(count)
reference = None

for index, subpath in sampling(paths):
    if not arguments.quiet:
        sys.stdout.write('\r[%d/%d]' % (index, count - 1))
        sys.stdout.flush()

    path = os.path.join(directory, subpath)
    model = biophysics.Model(path)
    error[index] = max(model.bond_lengths())

print

plot(models, error)

reference
plot((models[0], models[-1]), (error[0], error[0]), 'k:')

xlabel("Iteration")
ylabel(u"Longest Bond Length (Å)")

#xticks((0, 5000, 10000, 15000, 20000, 25000, 30000),
        #('0', '50K', '10K', '15K', '20K', '25K', '30K'))

if arguments.solutions:
    with open('solutions.dat') as file:
        solutions = file.readlines()
        solutions = map(int, solutions)

    plot(models, solutions)

if arguments.output:
    title("")
    savefig(arguments.output, dpi=arguments.dpi)

if not arguments.quiet and not os.fork(): 
    title(arguments.job.replace('_', '-'))
    show()
