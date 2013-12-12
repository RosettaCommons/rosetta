#!/usr/bin/env python

import os, re
import argparse
import numpy as np
import pylab as plt
from graphics import tango

parser = argparse.ArgumentParser()
parser.add_argument('jobs', nargs='+')
parser.add_argument('--include', '-i', default='.')
parser.add_argument('--exclude', '-x', default='^$')
#parser.add_argument('--limit-calc', '-lc', type=int, default=1000)
parser.add_argument('--limit-display', '-ld', type=int, default=5000)
parser.add_argument('--force', '-f', action='store_true')

arguments = parser.parse_args()

colors = {
    'uniform': tango.blue[1],
    'rama': tango.green[2],
    'biased-rama': tango.red[1]
}

for job in arguments.jobs:
    title = os.path.basename(job)
    trajectory_path = os.path.join(job, 'histogram/torsions.npz')
    correlation_path = os.path.join(job, 'histogram/correlations.npz')
    trajectories = np.load(trajectory_path)

    if os.path.exists(correlation_path) and not arguments.force:
        correlations = np.load(correlation_path)
    else:
        correlations = {}

    for name, trajectory in trajectories.items():
        if re.search(arguments.exclude, name): continue
        if not re.search(arguments.include, name): continue

        #trajectory = trajectory[:arguments.limit]

        print 'Thinking about %s/%s...' % (title, name)
        print '  %d steps' % len(trajectory)

        if name in correlations:
            correlation = correlations[name]
        else:
            correlation = np.correlate(trajectory, trajectory, mode='full')
            correlation = correlation[correlation.size/2:] / correlation.max()
            correlations[name] = correlation

        steps = np.arange(arguments.limit_display)
        correlation = correlation[:arguments.limit_display]
        #color = colors[title]
        label = '%s: %s' % (title, name)
        style = '--' if name.startswith('phi') else '-'

        plt.plot(steps, correlation, style, label=label) #, color=color)

    np.savez(correlation_path, **correlations)

plt.title("Autocorrelation")
plt.ylabel("Correlation")
plt.xlabel("Steps")
plt.legend(loc='upper right')
plt.ylim(-0.03, 1.0)

if not os.fork(): plt.show()

