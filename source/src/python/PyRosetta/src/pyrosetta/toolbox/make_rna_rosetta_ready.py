#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   make_rna_rosetta_ready.py
## @brief
## @author ?

from __future__ import print_function

import string
from sys import argv,stderr
from os import popen,system
from os.path import exists,dirname,basename,abspath

#assert( len(argv)>2)
pdbname = argv[1]

if (pdbname[-4:] != '.pdb' and pdbname[-7:] != '.pdb.gz'):
    pdbname += '.pdb'

removechain = 0
if argv.count('-nochain'):
    pos = argv.index('-nochain')
    del( argv[ pos ] )
    removechain = 1

ignore_chain = 0
if argv.count('-ignorechain'):
    pos = argv.index('-ignorechain')
    del( argv[ pos ] )
    ignore_chain = 1

chainids = []
if len( argv ) > 2:
    chainids = argv[2:]

if len(chainids) > 0 and len(chainids[0])==1:
    pdbnames = [ pdbname ]
else:
    pdbnames = argv[1:]
    ignore_chain = 1

fastaid = stderr
num_model = 0

max_model = 0 # for virus

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u'
              }

for pdbname in pdbnames:
    #netpdbname = '/net/pdb/' + pdbname[1:3] + '/' + pdbname
    #if not exists(netpdbname):
    netpdbname = pdbname

    print( 'Reading ... {}'.format(netpdbname) )
    if not exists( netpdbname ):
        print( 'DOES NOT EXIST: {}'.format(netpdbname) )
        continue

    if ( netpdbname[-3:] == '.gz' ):
        lines = popen( 'gzcat '+netpdbname ).readlines()
    else:
        lines = open(netpdbname,'r').readlines()

    outfile = string.lower( basename( pdbname ) )
    #outfile = dirname( abspath(outfile) ) + '/' + string.lower( basename(outfile) )
    outfile = outfile.replace( '.pdb', '_RNA.pdb').replace('.gz','');
    outid = open( outfile, 'w')

    print( 'Writing ... {}'.format(outfile) )

    #fastafile = pdbname[0:4]+chainid+'.pdb.fasta'
    #fastaid = open( fastafile, 'w')
    #print 'Writing ... '+fastafile
    fastaid.write('>'+pdbname+'\n');

    oldresnum = '   '
    count = 0;

    for i in range( len( chainids ) ) :
        if chainids[i] == '_':
            chainids[i] = ' '

    goodnames = [' rA',' rC',' rG',' rU']
    for line in lines:
        if len(line)>5 and line[:6]=='ENDMDL':
            num_model += 1
            if num_model > max_model:  break #Its an NMR model.
        if len(line) <= 21:  continue
        if (line[21] in chainids or ignore_chain):
            line_edit = line

            if line[0:3] == 'TER':
                continue
            elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
                if (line_edit[12:14] == 'SE'):
                    line_edit = line_edit[0:12]+' S'+line_edit[14:]
                if len(line_edit)>75:
                    if (line_edit[76:78] == 'SE'):
                        line_edit = line_edit[0:76]+' S'+line_edit[78:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='5BU'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='OMC'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='5MC'):
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]==' DC'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='CBR'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='CB2'):
                line_edit = 'ATOM  '+line[6:17]+'  C'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='2MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='H2U'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='PSU'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='5MU'):
                line_edit = 'ATOM  '+line[6:17]+'  U'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='OMG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='7MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='1MG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]==' YG'):
                line_edit = 'ATOM  '+line[6:17]+'  G'+line[20:]
            elif (line[0:6] == 'HETATM') & (line[17:20]=='1MA'):
                line_edit = 'ATOM  '+line[6:17]+'  A'+line[20:]

            #Don't save alternative conformations.
            if line[16] == 'A':
                continue;

            if line_edit[0:4] == 'ATOM':
                resnum = line_edit[23:26]
                if not resnum == oldresnum: #  or line_edit[12:16] == ' P  ':
                    longname = line_edit[17:20]
                    if longname == '  G':
                        longname = ' rG'
                    elif longname == '  A':
                        longname =   ' rA'
                    elif longname == '  C':
                        longname =   ' rC'
                    elif longname == '  U':
                        longname =   ' rU'
                    elif longname == 'G  ':
                        longname =   ' rG'
                    elif longname == 'A  ':
                        longname =   ' rA'
                    elif longname == 'C  ':
                        longname =   ' rC'
                    elif longname == 'U  ':
                        longname =   ' rU'
                    #elif longname == ' DG':
                    #    longname = ' rG'
                    #elif longname == ' DA':
                    #    longname = ' rA'
                    #elif longname == ' DC':
                    #    longname = ' rC'
                    #elif longname == ' DT':
                    #    longname = ' rU'
                    elif longname == 'GUA':
                        longname = ' rG'
                    elif longname == 'ADE':
                        longname = ' rA'
                    elif longname == 'CYT':
                        longname = ' rC'
                    elif longname == 'URA':
                        longname = ' rU'
                    elif longname == 'URI':
                        longname = ' rU'
                    else:
                        if longname not in goodnames:    continue

                    if longer_names.has_key(longname):
                        fastaid.write( longer_names[longname] );
                    else:
                        fastaid.write( 'X')

                    #print "AAH ==> " ,  resnum, oldresnum, line_edit
                    count = count + 1

                oldresnum = resnum

                if not longname in goodnames:
                    continue

                newnum = '%4d' % count
                line_edit = line_edit[0:16] + ' ' + longname + line_edit[20:22] + \
                            newnum + line_edit[26:]
                if removechain:
                    line_edit = line_edit[0:21]+'  '+line_edit[23:]

                line_edit = line_edit.replace( 'HO2\'', '2HO*' )
                line_edit = line_edit.replace( 'HO5\'', '5HO*' )
                line_edit = line_edit.replace( 'H5\'\'', '2H5*' )

                line_edit = line_edit.replace('\'','*')
                line_edit = line_edit.replace('OP1','O1P')
                line_edit = line_edit.replace('OP2','O2P')

                outid.write(line_edit)


    fastaid.write('\n')
    outid.close()
    #fastaid.close()
