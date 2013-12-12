#!/usr/bin/env python

import argparse
from os import fork
from sys import exit
from pylab import *

parser = argparse.ArgumentParser()
parser.add_argument('--iterations', '-i', type=int, default=10000)
parser.add_argument('--mean', '-c', type=float, default=0)
parser.add_argument('--stddev', '-s', type=float, default=1)
parser.add_argument('--no-counter-bias', '-ncb', action='store_true')
parser.add_argument('--no-score-function', '-nsf', action='store_true')
parser.add_argument('--plot-trajectory', '-pt', action='store_true')
parser.add_argument('--plot-histogram', '-ph', action='store_true')
parser.add_argument('--plot-correlation', '-pc', action='store_true')
arguments = parser.parse_args()

def propose_move(input):
    return normal(arguments.mean, arguments.stddev)

@vectorize
def proposal_prob(move):
    x = move - arguments.mean
    d = 2 * arguments.stddev**2
    return exp(-x**2 / d) / sqrt(pi * d)

def initial_accept(input, trial):
    if arguments.no_counter_bias:
        return True

    input_prob = proposal_prob(input)
    trial_prob = proposal_prob(trial)

    chance = min(1, input_prob / trial_prob)
    return chance > uniform()

def final_accept(input, trial):
    if arguments.no_score_function:
        return True

    input_prob = proposal_prob(input)
    trial_prob = proposal_prob(trial)

    chance = min(1, trial_prob / input_prob)
    return chance > uniform()


def plot_trajectory(trajectory):
    figure()
    title("Trajectory")
    ylabel("Move")
    xlabel("Steps")

    x = arange(len(trajectory))
    plot(x, trajectory)

def plot_histogram(trajectory):
    figure()
    title("Histogram")
    ylabel("Counts")
    xlabel("Move")

    x = linspace(-10, 10, num=500)
    plot(x, proposal_prob(x))

    bins = linspace(-10, 10, num=100)
    hist(trajectory, normed=True, bins=bins)

def plot_correlation(trajectory):
    figure()
    title("Autocorrelation")
    ylabel("Correlation")
    xlabel("Steps")

    x = arange(len(trajectory))
    result = correlate(trajectory, trajectory, mode='full')
    plot(x, result[result.size/2:] / max(result))


# Run the Monte Carlo simulation.

input = arguments.mean
trajectory = []

for i in range(arguments.iterations):
    trial = propose_move(input)

    if initial_accept(input, trial):
        if final_accept(input, trial):
            input = trial
    
    trajectory.append(input)

# Plot the resulting trajectory.

if arguments.plot_trajectory: plot_trajectory(trajectory)
if arguments.plot_histogram: plot_histogram(trajectory)
if arguments.plot_correlation: plot_correlation(trajectory)

if not fork(): show()
