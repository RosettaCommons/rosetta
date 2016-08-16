#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author: Rhiju Das

import string
from sys import argv,stdout
from os import popen,system
from os.path import exists,basename
from amino_acids import longer_names


pdbnames = argv[1:]

#chainid = ' '
#if len(argv)>2:
#    chainid = argv[2]

for pdbname in pdbnames:
#    if (pdbname[-4:] != '.pdb'):
#        pdbname += '.pdb'

    outfile = pdbname

    removechain = 0
    if argv.count('-nochain'):
        removechain = 1

    netpdbname = pdbname
    assert( exists(netpdbname))
    #print 'Reading ... '+netpdbname

    lines = open(netpdbname,'r').readlines()

    #outid = open( outfile, 'w')
    #print 'Writing ... '+pdbname

    #fastafile = pdbname+'.fasta'
    #fastaid = open( fastafile, 'w')
    #print 'Writing ... '+fastafile

    fastaid = stdout

    fastaid.write('>'+basename(pdbname)+'\n');

    oldresnum = '   '
    count = 0;
    for line in lines:
        if (len(line)>20): # and (chainid == line[21]):
            line_edit = line
            if line[0:3] == 'TER':
                fastaid.write('\n')
                #    break
            elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
                if (line_edit[12:14] == 'SE'):
                    line_edit = line_edit[0:12]+' S'+line_edit[14:]
                if len(line_edit)>75:
                    if (line_edit[76:78] == 'SE'):
                        line_edit = line_edit[0:76]+' S'+line_edit[78:]

            if line_edit[0:4] == 'ATOM':
                resnum = line_edit[23:26]
                if not resnum == oldresnum:
                    count = count + 1
                    longname = line_edit[17:20]
                    if longer_names.has_key(longname):
                        fastaid.write( longer_names[longname] );
                    else:
                        fastaid.write( 'X')
                oldresnum = resnum

                newnum = '%3d' % count
                line_edit = line_edit[0:23] + newnum + line_edit[26:]
                if removechain:
                    line_edit = line_edit[0:21]+' '+line_edit[22:]

                #outid.write(line_edit)
    fastaid.write('\n')


    #outid.close()
    #fastaid.close()
