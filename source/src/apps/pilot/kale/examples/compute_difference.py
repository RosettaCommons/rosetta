#!/usr/bin/env python

differences = [[], []]
group = -1

with open('comparison.dat') as file:
    for line in file:
        try:
            a, b = line.strip().split()
            a = float(a)
            b = float(b)

            #print '%7.3f   %8.4f        %8.4f' % (a, b, a - b)
            differences[group].append(a - b)

        except ValueError:
            #print line.strip()

            if line.startswith('='):
                group += 1

from os import fork
from pylab import *

hist(differences[0], label='Solution 5', alpha=0.5, bins=20)
hist(differences[1], label='Solution 6', alpha=0.5, bins=20)

if not fork(): show()
