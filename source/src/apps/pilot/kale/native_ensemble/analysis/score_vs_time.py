#!/usr/bin/env python

import argparse
import helpers, schema
import numpy, pylab
from cStringIO import StringIO

parser = helpers.define_args('plot')
arguments = parser.parse_args()
job = arguments.job

with helpers.connect_to_database(arguments.database):
    iterations = helpers.load_array(job, 'iterations')[1:]
    scores = helpers.load_array(job, 'scores')[1:]

pylab.title('%s' % helpers.job_title(arguments))
pylab.plot(iterations, scores)
pylab.xlabel('Iterations')
pylab.ylabel('Score')

helpers.make_plot(arguments)


