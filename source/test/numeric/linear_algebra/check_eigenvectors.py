#!/usr/bin/env python

import sys
from scipy import *
from scipy.linalg import *

print "Paste the output from show_matrices(), and press Ctrl-D to analyze..."
exec sys.stdin
print; print 79 * '-'; print

a = matrix(A)
b = matrix(B)
x = array(x)
x = x[...,0] + 1j * x[...,1]

print "The columns of this matrix should be the eigenvalues..."
print array(a * x) / array(b * x)
print

print "The following numbers should all be unity..."
for col in range(x.shape[0]):
    print norm(x[:,col])
