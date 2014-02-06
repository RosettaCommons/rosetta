#!/usr/bin/env python

import helpers
import numpy, pylab

parser = helpers.define_args('torsion', 'plot')
parser.add_argument('--ymin', '-y', type=float, default=0)
parser.add_argument('--window', '-w', type=float)
arguments = parser.parse_args()

def find_nonmodulo_autocorrelation(x):
    """ The standard correlation function isn't meant for modulo variables, 
    like dihedral angles.  I'm sure this version is much faster, but it will 
    choke on angles constantly jumping between 180 and -180. """
    x -= numpy.mean(x)
    x /= numpy.linalg.norm(x)
    return numpy.correlate(x, x, mode='full')[len(x)-1:]

def find_autocorrelation(x):
    """ This homebrew autocorrelation function is specifically written to 
    handle modulo variable like dihedral angles.  Think of the cosine function 
    as an alternative, and more appropriate, metric for overlap. """
    from numpy import ones, sum, cos, radians
    n = len(x); c = ones(n)
    for dt in range(1, n):
        c[dt] = sum(cos(radians(x[:-dt] - x[dt:]))) / (n - dt)
    return c


with helpers.connect_to_database(arguments.database):
    iterations = helpers.load_array(arguments.job, 'iterations')
    torsions = helpers.load_arrays(arguments.job, 'torsions')

for key in helpers.pick_torsions(arguments, torsions):
    autocorrelation = find_autocorrelation(torsions[key])
    color = pylab.cm.Dark2(numpy.median(autocorrelation))
    pylab.plot(iterations, autocorrelation, label=key, color=color)

pylab.title('Autocorrelation: Job %d' % arguments.job)
pylab.xlabel('iterations')
pylab.ylabel('correlation')
pylab.ylim(arguments.ymin, 1)
if arguments.window: pylab.xlim(0, arguments.window)

helpers.make_plot(arguments)
