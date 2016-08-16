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

#score output as defined by the tags in the first line
resultfile = './output/enz_score.out'

isTestPassed = True

lines = open( resultfile ).readlines()

num_lines=0
tot_recovery=0
target_recovery=0.38

#gets the tags (first line in the output) and figures out
#column position for the sequence recovery
cols0=lines[0].split()
index_for_tot_seq_rec= cols0.index("tot_seq_recovery")

#gets specific seq. recovery for each protein and ads them up
for line in lines[1:]:
    try:
        cols = line.split()

        recovery = float( cols[index_for_tot_seq_rec] )
        #write individual seq. recoveries in the log
        output = ' %6.2f ' % ( recovery )
        fid.write(output)
        tot_recovery = tot_recovery+recovery
        num_lines=num_lines+1
    except ValueError: pass

tot_recovery=tot_recovery/num_lines

##print tot_recovery
my_output = 'seq. recovery: %6.2f   [should be more than %6.2f]' % ( tot_recovery, target_recovery )

fid.write( my_output  )
fid.close()

isTestPassed &= tot_recovery > target_recovery
fid_yaml.write( "'%s' : %6.2f, " % ( "seq_recovery", tot_recovery ) )
fid_yaml.write( "'_isTestPassed' : %s, " % isTestPassed )
fid_yaml.write(' }\n')
fid_yaml.close()



