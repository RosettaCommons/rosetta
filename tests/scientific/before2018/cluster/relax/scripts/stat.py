#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import string
import sys
import math

f = sys.stdin.read()
lines = string.split(f,'\n')

data = []
sum = 0.0
sumcnt = 0
sqrsum = 0
n = 0.0
for line in lines:
  token = string.split(line)
  if len(token)<=0: continue
  d = float( token[0] )
  data.append( d )
  sum += d
  sqrsum += d * d
  n = n + 1.0

if( n > 0): var = sqrsum / n - sum * sum / ( n* n)
else:       var = 0;

if( var >= 0 ): stdev = math.sqrt(var)
else:          stdev = 0;

data.sort()

#print "sum sqrsum av var stdev  0.001 0.01 0.1   0.5   0.9 0.99 0.999"
print sum / n, var, stdev, int(n) ,"   ", sum, sqrsum, "   ", data[ int(n * 0.001) ] ,  \
      data[ int(n * 0.01) ], data[ int(n*0.1) ] , "   ", data[ int(n * 0.5) ] ,  "   ",  \
			data[ int(n*0.9) ], data [ int(n*0.99) ], data [ int(n*0.999) ]
