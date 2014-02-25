#!/usr/bin/env python
# encoding: utf-8

import os
import argparse
import itertools

import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as fancy_axes
import mpl_toolkits.axisartist as axis_artist

from pylab import *
from numpy import *
from pprint import pprint
from graphics import tango

parser = argparse.ArgumentParser()
parser.add_argument('input')
arguments = parser.parse_args()

path = os.path.join(arguments.input)

pivots_to_keys = {}
keys_to_pivots = {}

iterator = itertools.combinations(range(1, 13), 2)
combos = [tuple(x) for x in iterator if x[0] + 1 != x[1]]

length_key = lambda x: x[1] - x[0]
combos.sort(key=length_key)
lengths = map(length_key, combos)

for key, pivots in sampling(combos):
    pivots_to_keys[pivots] = key
    keys_to_pivots[key] = pivots

num_combos = len(pivots_to_keys)
histograms = [[] for x in range(num_combos)]

def pivots_from_line(line):
    fields = line.split()
    return int(fields[2]), int(fields[3])

def key_from_pivots(pivots):
    return pivots_to_keys[pivots]


with open(path) as file:
    for line in file:
        if line.startswith('PIVOTS CHOOSING'):
            pivots = pivots_from_line(line)
            key = key_from_pivots(pivots)
            histograms[key].append(0)

        if line.startswith('PIVOTS ITERATING'):
            pivots = pivots_from_line(line)
            key = key_from_pivots(pivots)
            histograms[key][-1] += 1

axes = subplot(111)
twin = axes.twiny()

top_labels, top_label_locations = [], []
bottom_labels, bottom_label_locations = [], []

for i in range(lengths[-1] + 1):
    try: x = lengths.index(i)
    except ValueError: continue
    bottom_labels.append(u'%s' % (i+1))
    bottom_label_locations.append(x+1)

for key, pivots in keys_to_pivots.items():
    top_labels.append('%d-%d' % pivots)
    top_label_locations.append(key + 1)

artists = axes.boxplot(histograms)
black = tango.black

plt.setp(artists['boxes'], color=black)
plt.setp(artists['medians'], color=black, linewidth=3, solid_capstyle='butt')
plt.setp(artists['whiskers'], color=black, linestyle='-')
plt.setp(artists['caps'], color=black)
plt.setp(artists['fliers'], color=black, marker='', markersize=1)

for box in artists['boxes']:
    coordinates = zip(box.get_xdata(), box.get_ydata())
    polygon = plt.Polygon(coordinates, facecolor=tango.brown[0])
    axes.add_patch(polygon)

xpadding = 0.02 * num_combos
axes.set_xlim(1 - xpadding, num_combos + xpadding)
axes.set_xticks(bottom_label_locations)
axes.set_xticklabels(bottom_labels, size='small')
axes.set_xlabel('Residues Between Pivots')

twin.set_xticks(top_label_locations)
twin.set_xticklabels(top_labels, rotation='vertical', size='x-small')
twin.set_xlim(*axes.get_xlim())

ymin, ymax = axes.get_ylim()
ypadding = 0.03 * ymax
axes.set_ylim(1 - ypadding, min(ymax, 2000) + ypadding)
axes.set_ylabel('Sampling Attempts')
plt.setp(axes.get_yticklabels(), size='small')

plt.savefig('subloop-histogram.png', dpi=150)
if not os.fork():
    plt.show()

