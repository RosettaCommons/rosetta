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
fid2 = open( 'staResult', 'w' )

fid_yaml = open( '.results.yaml', 'w')
fid_yaml.write('{ ')

resultfile = 'inputs/2r8s_RNA.sequence_recovery.txt'

titles = { 'SINGLE_STRANDED':'no base pairs', 'DOUBLE_STRANDED':'Watson/Crick bp', 'TERTIARY_STRUCT':'non-Watson/Crick bp','OVERALL':'all RNA'}

titles_yaml = { 'SINGLE_STRANDED':'RNA_SeqRecoveryNoBasePairs', 'DOUBLE_STRANDED':'RNA_SeqRecoveryWC_BasePairs', 'TERTIARY_STRUCT':'RNA_SeqRecoveryNonWC_BasePairs','OVERALL':'RNA_SeqRecoveryOverall'}

target_recoveries = {'SINGLE_STRANDED':0.30, 'DOUBLE_STRANDED':0.35, 'TERTIARY_STRUCT':0.7, 'OVERALL':0.45}

isTestPassed = True

lines = open( resultfile ).readlines()

for line in lines[-4:]:
    cols = line.split()

    tag = cols[0]
    recovery = float( cols[-1] )

    my_output = '\n%20s, seq. recovery: %6.2f   [should be more than %6.2f]' % ( titles[ tag ], recovery, target_recoveries[ tag ] )
    fid.write( my_output  )
    fid2.write( my_output )

    fid_yaml.write( "'%s' : %6.2f, " % ( titles_yaml[tag], recovery ) )

    isTestPassed &= recovery > target_recoveries[ tag ]


fid_yaml.write( "'_isTestPassed' : %s, " % isTestPassed )

fid.close()
fid2.close()

fid_yaml.write(' }\n')
fid_yaml.close()
