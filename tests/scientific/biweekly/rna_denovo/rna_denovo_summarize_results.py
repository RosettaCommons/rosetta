#!/usr/bin/python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

from os.path import exists
from os import getcwd
import string
from sys import exit,stdout

fid = open( '.results.log', 'w' )

fid_yaml = open( '.results.yaml', 'w')
fid_yaml.write('{ ')

styles = ['FARNA','Hi-res refine' ]
files = ['outputs/chunk002_1lnt_LORES.out', 'outputs/chunk002_1lnt.out']
target_rms = [3.0, 2.75]

isTestPassed = True

for i in range( len( files )):
    file = files[ i ]

    if not exists( file ):
        fid.write( 'Hey, where is ', file, '?\n' )
        continue

    lines = open( file ).readlines()

    line = lines[1]
    cols = line.split()
    rms_index = cols.index( 'rms' )

    scores = []
    for line in lines[2:]:
        cols = line.split()
        if len( cols ) > 0 and cols[0] == 'SCORE:':
            try:
                scores.append( float( cols[ rms_index] ) )
            except:
                continue

    numscores = len( scores )
    scores.sort()
    median_score = scores[ numscores/4 ]
    fid.write( '\n%20s, 25th percentile rms (%d decoys): %6.2f   [should be less than %6.2f]' % ( styles[i], numscores, median_score, target_rms[i] ) )
    fid_yaml.write( "'%s_25thPercentileRMS' : %6.2f, " % ( styles[i], median_score ) )

    isTestPassed &= median_score < target_rms[i]


fid_yaml.write( "'_isTestPassed' : %s, " % isTestPassed )

fid.close()
fid_yaml.write(' }\n')
fid_yaml.close()
