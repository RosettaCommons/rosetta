#!/usr/bin/env python

import helpers
import numpy, pylab
from scipy.stats import gaussian_kde

parser = helpers.define_args('plot')
arguments = parser.parse_args()
job = arguments.job

with helpers.connect_to_database(arguments.database):
    rmsds = helpers.load_array(job, 'rmsds')[1:]
    scores = helpers.load_array(job, 'scores')[1:]

# If there aren't so many frames, just use one color.  Otherwise, calculate a 
# color-map by making a smoothed version of the 2D histogram.  The particular 
# method used here is a gaussian kernel density estimation.  Note that this 
# method tends to over-smooth multi-modal distributions.

frames = len(scores)

if frames < 10000:
    colors = numpy.zeros(frames)
else:
    skip = frames // 500
    estimator = gaussian_kde([rmsds[::skip], scores[::skip]])
    density = estimator.evaluate([rmsds, scores])
    colors = 2**density.clip(numpy.median(density), numpy.Inf)

style = {}
style['s'] = 10         # size
style['c'] = colors     # color
style['edgecolors'] = 'none'
style['cmap'] = pylab.cm.jet

pylab.title('%s' % helpers.job_title(arguments))
pylab.scatter(rmsds, scores, **style)
pylab.axvline(1.0, color='k', linestyle='--', zorder=0)
pylab.xlabel('RMSD')
pylab.ylabel('Score')
pylab.xlim(0, 6)

helpers.make_plot(arguments)
