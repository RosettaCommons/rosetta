#!/usr/bin/env python
# encoding: utf-8

import os, glob
import argparse
from pylab import *

parser = argparse.ArgumentParser()
parser.add_argument('jobs', nargs='+')
parser.add_argument('--top', '-t', type=int, default=500)
parser.add_argument('--source', '-s', default='old')
parser.add_argument('--ylim', nargs=2, type=int)
arguments = parser.parse_args()

def plot_job(job, color):
    launcher = os.path.join(job, 'job.sh')
    cluster_jobs = os.path.join(job, '*')

    if os.path.exists(launcher):
        subjobs = [job]
    else:
        subjobs = glob.glob(cluster_jobs)

    rmsds = []
    scores = []

    for subjob in subjobs:
        if arguments.source == 'old':
            parse_old_source(subjob, rmsds, scores)
        else:
            parse_new_source(subjob, rmsds, scores)

    if not scores or not rmsds:
        raise SystemExit

    size = 4 if len(rmsds) < 1000 else 1

    scatter(rmsds, scores,
            s=size, color=color, edgecolors='none', zorder=2, label=job)

def parse_new_source(subjob, rmsds, scores):
    path = os.path.join(subjob, 'score_vs_rmsd.txt')

    try:
        with open(path) as file:
            lines = file.readlines()
            lines = lines[:arguments.top]

        for line in lines:
            rmsd, score = map(float, line.split())
            rmsds.append(rmsd)
            scores.append(score)

    except IOError:
        print 'No data for:', path

def parse_old_source(subjob, rmsds, scores):
    path = os.path.join(subjob, 'output.txt')

    rmsd_header = 'protocols.loop_build.LoopBuildMover: loop_rms'
    score_header = 'protocols.loop_build.LoopBuildMover: total_energy'

    try:
        with open(path) as file:
            file.seek(0, 2)             # Go to the end of the file.
            file_size = file.tell()     # Figure out how big the file is.

            tail_start = max(file_size - 1000, 0)
            file.seek(tail_start, 0)
            lines = file.readlines()

        for line in lines:
            if line.startswith(rmsd_header):
                fields = line.split()
                rmsds.append(float(fields[-1]))
            if line.startswith(score_header):
                fields = line.split()
                scores.append(float(fields[-1]))

    except IOError:
        print 'No data for:', path


colors = [
        '#204a87',  # Blue
        '#a40000',  # Red
        '#4e9a06',  # Green
        '#ce5c00',  # Orange
        '#5c3566']  # Purple

for job, color in zip(reversed(arguments.jobs), colors):
    plot_job(job, color)

axvline(1, linestyle='--', color='gray', zorder=1)
legend()

xlabel(u'Backbone RMSD (Å)')
ylabel('Score (REU)')
xlim(0, 6)
if arguments.ylim: ylim(*arguments.ylim)

if not os.fork():
    show()
