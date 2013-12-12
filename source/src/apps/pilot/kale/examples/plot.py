#!/usr/bin/env python

import argparse
from graphics import tango
from pylab import *
from os import fork

parser = argparse.ArgumentParser()
parser.add_argument('--distributions', '-d', nargs='+', default=[])
parser.add_argument('--functions', '-f', nargs='+', default=[])
parser.add_argument('--scores', '-s', nargs='+', default=[])
parser.add_argument('--output', '-o')
arguments = parser.parse_args()

for path in arguments.functions + arguments.scores:     # (fold)
    data = {}

    # Read the function in from a text file.

    with open(path) as file:
        for line in file:
            phi, psi, value = line.split()
            phi, psi, value = int(phi), int(psi), float(value)

            data.setdefault(phi, {})[psi] = value

    # Convert the function into a 2D numpy array.

    phis = sort(data.keys())
    psis = sort(data.values()[0].keys())

    num_phi = len(phis)
    num_psi = len(psis)

    function = zeros((num_phi, num_psi))
    for i, phi in enumerate(phis):
        for j, psi in enumerate(psis):
            function[i,j] = data[phi][psi]

    # Normalize and plot the function.

    if path in arguments.scores:
        function = exp(-function)

    #imshow(function, interpolation='nearest')

    phi_function = function.sum(1) / function.sum()
    psi_function = function.sum(0) / function.sum()

    plot(phis, phi_function, label=path+' $(\phi)$', color=tango.red[1]) 
    plot(psis, psi_function, label=path+' $(\psi)$', color=tango.blue[1]) 

for path in arguments.distributions: # (fold)
    phi_hist = []
    psi_hist = []

    with open(path) as file:
        for line in file:
            phi, psi = line.split()
            phi_hist.append(float(phi))
            psi_hist.append(float(psi))

    #phi_counts, phi_bins = histogram(phi_hist, bins=360, normed=True)
    #psi_counts, psi_bins = histogram(psi_hist, bins=360, normed=True)
    bins = arange(360)
    function, phi_bins, psi_bins = histogram2d(phi_hist, psi_hist, bins)
    function /= function.max()

    phi_bins = phi_bins[:-1] + (phi_bins[1] - phi_bins[0])
    psi_bins = psi_bins[:-1] + (psi_bins[1] - psi_bins[0])

    imshow(function, interpolation='nearest')

    #plot(phi_bins, phi_counts, ':', color=tango.red[1])
    #plot(psi_bins, psi_counts, ':', color=tango.blue[1])


xlim(-180, 180)
title("Mismatch Between Rama Scoring & Rama Sampling")
ylabel("Phi")
xlabel("Psi")
legend(loc='upper left')
#colorbar()

if arguments.output:
    savefig(arguments.output)
elif not fork():
    show()





