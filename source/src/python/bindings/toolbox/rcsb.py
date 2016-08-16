#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   rcsb.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

import os
import urllib

from rosetta import Pose
from rosetta import pose_from_file

# other tools
from cleaning import cleanATOM
from cleaning import cleanCRYS

# retreives pdbData from rcsb for  <pdb_code>
# ADD NAMING OPTION
def load_from_rcsb( pdb_code , pdb_outfile = '' ):
    """
    Writes PDB data for RCSB data for  <pdb_code>  into the file  <pdb_code>.pdb

    example:
        load_from_rcsb('1YY8')
    See also:
        Pose
        pose_from_file
        pose_from_rcsb
        pose_from_sequence
        cleanATOM
        cleanCRYS
    """
    if pdb_code:    # if something input...
        pdb_code = pdb_code.upper()
        try:
            filename = urllib.urlretrieve('http://www.rcsb.org/pdb/files/' + pdb_code + '.pdb')[0]
        except:
            raise IOError('Cannot access the PDB database, please check your Internet access')
        else:
            if (os.path.getsize(filename) > 1500):    # arbitrary 1500...then pdb_code was invalid
                # load in the data
                pdb_file = open(filename)
                pdb_data = pdb_file.readlines()
                pdb_file.close()

                # setup proper naming
                pdb_code = pdb_code + '.pdb'
                if not pdb_outfile:
                    # default the name to the <pdb_code>.pdb
                    pdb_outfile = pdb_code
                if os.path.exists( os.getcwd() + '/' + pdb_outfile ):
                    print 'the file',pdb_outfile,'already exists, this file will be overwritten'
                #if input('Do you want to overwrite ' + pdbCode + '.pdb')
                pdb_file = open(pdb_outfile,'w')
                pdb_file.writelines(pdb_data)
                pdb_file.close()

                print 'PDB',pdb_code[:-4],'successfully loaded from rcsb into',pdb_outfile
#                if auto_clean:
#                    cleanATOM(pdb_code)
            else:
                raise IOError('Invalid PDB code')
        os.remove(filename)    # remove temp file

# packaged my method to fit naming
def pose_from_rcsb( pdb_code , ATOM = True , CRYS = False , pdb_outfile = '' ):
    """
    Returns a pose for RCSB PDB  <pdb_code> , also writes this data to
    <pdb_code>.pdb, optionally calls cleanATOM and cleanCYRS

    example:
        pose=pose_from_rcsb('1YY8')
    See also:
        Pose
        pose_from_file
        pose_from_sequence
        load_from_rcsb
        cleanATOM
        cleanCRYS
    """
    load_from_rcsb(pdb_code,pdb_outfile)
    # ensure the names are proper
    edit = -4
    if not pdb_outfile:
        pdb_outfile = pdb_code + '.pdb'
    else:
        edit = 0
#    if not pdb_outfile[:-4]=='.pdb':
#        pdb_outfile = pdb_outfile
    # cleaning calls
    if ATOM:
        cleanATOM(pdb_outfile,edit=edit)
        pdb_outfile = pdb_outfile[:edit]+'.clean.pdb'
    if CRYS:
        cleanCRYS(pdb_outfile)
        pdb_outfile = pdb_outfile[:edit]+'.mono.pdb'
    pose = pose_from_file(pdb_outfile)
    return pose

# retreives pdbData from rcsb for  <pdb_code>
# ADD NAMING OPTION
def load_fasta_from_rcsb( pdb_code , fasta_outfile ):
    if pdb_code:    # if something input...
        pdb_code = pdb_code.upper()
        try:
            filename = urllib.urlretrieve('http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=' + pdb_code)[0]
        except:
            raise IOError('Cannot access the PDB database, please check your Internet access')
        else:
            if (os.path.getsize(filename)):    # arbitrary 1500...then pdb_code was invalid
                pdb_file = open(filename)
                pdb_data = pdb_file.readlines()
                pdb_file.close()

#                pdb_code = pdb_code + '.fa'
                if not fasta_outfile:
                    fasta_outfile = pdb_code + '.fa'
                if os.path.exists( os.getcwd() + '/' + fasta_outfile ):
                    print 'the file',fasta_outfile,'already exists, this file will be overwritten'
                #if input('Do you want to overwrite ' + pdbCode + '.pdb')
                pdb_file = open(fasta_outfile,'w')
                pdb_file.writelines(pdb_data)
                pdb_file.close()

                print 'PDB',pdb_code,'sequence successfully loaded from rcsb into',fasta_outfile
            else:
                raise IOError('Invalid PDB code')
        os.remove(filename)    # remove temp file


