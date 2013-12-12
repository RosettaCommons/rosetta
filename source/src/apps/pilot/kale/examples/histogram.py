#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bins', '-b', type=int, default=360//5)
arguments = parser.parse_args()


def plot_rama_kic():
    phis = []
    psis = []

    with open('rama.dat') as file:
        for line in file:
            phi, psi = line.split()
            phis.append(float(phi) % 360)
            psis.append(float(psi) % 360)

    bins = np.linspace(0, 360, num=arguments.bins)
    axis = bins[1:] - (bins[1] - bins[0]) / 2

    counts = np.histogram(phis, bins=bins, density=True)[0]
    plt.plot(axis, counts, label='phi rama')

    counts = np.histogram(psis, bins=bins, density=True)[0]
    plt.plot(axis, counts, label='psi rama')

def plot_uniform_kic():
    phis = []
    psis = []

    with open('uniform.dat') as file:
        for line in file:
            phi, psi = line.split()
            phis.append(float(phi) % 360)
            psis.append(float(psi) % 360)

    bins = np.linspace(0, 360, num=arguments.bins)
    axis = bins[1:] - (bins[1] - bins[0]) / 2

    counts = np.histogram(phis, bins=bins, density=True)[0]
    plt.plot(axis, counts, label='phi uniform')

    counts = np.histogram(psis, bins=bins, density=True)[0]
    plt.plot(axis, counts, label='psi uniform')


def plot_ideal_torsions():
    dist = np.zeros((360, 360))

    with open('rama/score/correct.dat') as file:
        for line in file:
            phi, psi, prob = line.split()
            phi = int(phi) % 360
            psi = int(psi) % 360
            dist[phi,psi] = float(prob)

    dist = np.exp(-dist)
    axis = np.arange(0, 360)

    phis = dist.sum(1) / dist.sum()
    psis = dist.sum(0) / dist.sum()

    plt.plot(axis, phis, label='phi expected')
    plt.plot(axis, psis, label='psi expected')


plot_rama_kic()
plot_uniform_kic()
plot_ideal_torsions()

plt.xlim(0, 360)
plt.legend()
plt.show()
