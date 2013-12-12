#!/usr/bin/env python

from math import *
from numpy import *
from numpy.linalg import *

import argparse
from biophysics import pdb

parser = argparse.ArgumentParser()
parser.add_argument('pdb', nargs=2)
arguments = parser.parse_args()
backbones = []

def length(model, index):
    a = model.coordinates[index % wrap]
    b = model.coordinates[(index + 1) % wrap]

    return norm(a - b)

def angle(model, index):
    a = model.coordinates[(index - 1) % wrap]
    b = model.coordinates[(index + 0) % wrap]
    c = model.coordinates[(index + 1) % wrap]

    r = (a - b) / norm(a - b)
    s = (c - b) / norm(c - b)

    return degrees(acos(dot(r, s)))

def torsion(model, index):
    a = model.coordinates[(index - 1) % wrap]
    b = model.coordinates[(index + 0) % wrap]
    c = model.coordinates[(index + 1) % wrap]
    d = model.coordinates[(index + 2) % wrap]

    q = b - a
    r = c - b
    s = d - c

    u = norm(r) * dot(q, cross(r, s))
    v = dot(cross(q, r), cross(r, s))

    return degrees(atan2(u, v)) % 360


for path in arguments.pdb:
    loop = pdb.Model(path)
    backbone = loop.select_atoms('N', 'CA', 'C')
    backbones.append(backbone)

wrap = len(backbones[0].coordinates)
indices = backbones[0].indices()

print "Lengths"
print "======="
for index in indices:
    lengths = [length(backbone, index) for backbone in backbones]
    readout = lengths[0], lengths[1], lengths[0] - lengths[1]
    print '%5.3f  %5.3f  %6.3f' % readout
print

print "Angles"
print "======"
for index in indices:
    angles = [angle(backbone, index) for backbone in backbones]
    readout = angles[0], angles[1], angles[0] - angles[1]
    print '%7.3f  %7.3f  %7.3f' % readout
print

print "Torsions"
print "========"
for index in indices:
    torsions = [torsion(backbone, index) for backbone in backbones]
    readout = torsions[0], torsions[1], torsions[0] - torsions[1]
    print '%7.3f  %7.3f  %7.3f' % readout
print
