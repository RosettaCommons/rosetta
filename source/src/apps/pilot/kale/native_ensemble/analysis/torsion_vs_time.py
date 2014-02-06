#!/usr/bin/env python

import helpers
import numpy, pylab

parser = helpers.define_args('torsion', 'plot')
arguments = parser.parse_args()

with helpers.connect_to_database(arguments.database):
    iterations = helpers.load_array(arguments.job, 'iterations')
    torsions = helpers.load_arrays(arguments.job, 'torsions')

for key in helpers.pick_torsions(arguments, torsions):
    pylab.plot(iterations, torsions[key] % 360, '.', label=key)

pylab.title('Torsion vs. Time: Job %d' % arguments.job)
pylab.xlabel('Iterations')
pylab.ylabel('Torsion (deg)')
pylab.ylim(0, 360)
pylab.yticks(range(0, 361, 60))

helpers.make_plot(arguments)
