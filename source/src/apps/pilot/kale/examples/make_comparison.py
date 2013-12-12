#!/usr/bin/env python

import re
import pprint
import argparse

from math import *
from numpy import *
from numpy.linalg import *
from biophysics import pdb

parser = argparse.ArgumentParser()
parser.add_argument('solution', type=int, choices=range(1,7))
arguments = parser.parse_args()

expected = {
        'torsions': [],
        'angles': [],
        'lengths': [] }

measured = {
        'torsions': [],
        'angles': [],
        'lengths': [] }

# Parse expected geometry

collecting_data = False
current_section = None
section_pattern = re.compile(r'(\w+) \((\d)/\d\)')

with open('expected.dat') as file:
    for line in file:
        try:
            value = float(line)
            if collecting_data:
                expected[current_section].append(value)

        except ValueError:
            match = section_pattern.match(line)
            if match:
                collecting_data = match.group(2) == str(arguments.solution)
                current_section = match.group(1).lower()

expected['torsions'] = array(expected['torsions'][1:-2])
expected['angles'] = array(expected['angles'][1:-2])
expected['lengths'] = array(expected['lengths'][1:-2])

# Parse measured geometry.

def length(index):
    a = coordinates[index % wrap]
    b = coordinates[(index + 1) % wrap]

    return norm(a - b)

def angle(index):
    a = coordinates[(index - 1) % wrap]
    b = coordinates[(index + 0) % wrap]
    c = coordinates[(index + 1) % wrap]

    r = (a - b) / norm(a - b)
    s = (c - b) / norm(c - b)

    return degrees(acos(dot(r, s)))

def torsion(index):
    a = coordinates[(index - 1) % wrap]
    b = coordinates[(index + 0) % wrap]
    c = coordinates[(index + 1) % wrap]
    d = coordinates[(index + 2) % wrap]

    q = b - a
    r = c - b
    s = d - c

    u = norm(r) * dot(q, cross(r, s))
    v = dot(cross(q, r), cross(r, s))

    return degrees(atan2(u, v)) % 360


loop = pdb.Model('solution.%d.pdb' % arguments.solution)
backbone = loop.select_atoms('N', 'CA', 'C')
coordinates = backbone.get_coordinates()
indices = backbone.indices()
wrap = len(coordinates)

measured['torsions'] = array([torsion(x) for x in indices][1:-5])
measured['angles'] = array([angle(x) for x in indices][1:-5])
measured['lengths'] = array([length(x) for x in indices][1:-5])

# Output the results.

torsion_error = measured['torsions'] - expected['torsions']
angle_error = measured['angles'] - expected['angles']
length_error = measured['lengths'] - expected['lengths']

print 'Torsion Error:', max(abs(torsion_error))
print 'Angle Error:  ', max(abs(angle_error))
print 'Length Error: ', max(abs(length_error))

