#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## Authors: Mike Tyka

import sys
import string

def Help():
    print sys.argv[0]+' pdbfile  tagfile'
    print
    sys.exit()


if len(sys.argv) < 3:
    Help()

pdbfile = sys.argv[1]
tagfile = sys.argv[2]

data = open(tagfile,'r').read()
lines = string.split(data,'\n')
tags = []
sumcolor = 0.0
for i in lines:
	token = string.split(i)
	if len(token) < 2 : continue;
	tags.append( [ int( float(token[0]) ),  float(token[1]) ] )
	sumcolor = sumcolor + float(token[1])

for i,j in enumerate(tags):
	tags[i][1] = j[1] / sumcolor * len(tags);


data = open(pdbfile,'r').read()
lines = string.split(data,'\n')
for i in lines:
	token = string.split(i)
	if len(token) < 6 : continue;
	if token[0] == "ATOM":
		res = int( i[ 6:11] )
		## find residue number
		temp = 0.0
		for j in tags:
		  if int(j[0]) == int(res):
				temp = j[1]

		temps = "% 4.2f"%float(temp);
		i = i[0:61] + temps + i[66:]

	print i


