#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   rcsb.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

from __future__ import print_function

import os, sys
import urllib

import pyrosetta.rosetta as rosetta

from pyrosetta import Pose, pose_from_file

# other tools
from pyrosetta.toolbox.cleaning import cleanATOM, cleanCRYS


if sys.version_info >= (3, 0):
    import urllib.request
    urllib_urlretrieve = urllib.request.urlretrieve

else:
    urllib_urlretrieve = urllib.urlretrieve


# retreives pdbData from rcsb for  <pdb_code>
# ADD NAMING OPTION
def load_from_rcsb(pdb_code, pdb_filename = None):
    """
    Writes PDB data for RCSB data for <pdb_code> into <pdb_filename>. If not
    specified, outputs file to <pdb_code>.pdb.

    Example:
        load_from_rcsb("1YY8")
    See also:
        Pose
        pose_from_file
        pose_from_rcsb
        pose_from_sequence
        cleanATOM
        cleanCRYS
    """
    pdb_code = pdb_code.upper()
    try:
        temp = urllib_urlretrieve("http://www.rcsb.org/pdb/files/" + pdb_code + ".pdb")[0]
    except:
        raise IOError("Cannot access the PDB database, please check your Internet access.")
    else:
        if (os.path.getsize(temp) > 1500):
            # Arbitrarily 1500... else pdb_code was invalid.
            with open(temp) as f:
                pdb_data = f.readlines()

            if not pdb_filename: pdb_filename = pdb_code + ".pdb"

            if os.path.exists(os.getcwd() + '/' + pdb_filename): print( "The file", pdb_filename, "already exists; this file will be overwritten." )

            with open(pdb_filename, 'w') as f: f.writelines(pdb_data)

            print( "PDB", pdb_code, "successfully loaded from the RCSB into", pdb_filename + '.' )
            #if auto_clean:
            #    cleanATOM(pdb_filename)
        else:
            raise IOError("Invalid PDB code")
        os.remove(temp)  # Remove temp file.


def pose_from_rcsb(pdb_code, ATOM = True, CRYS = False):
    """
    Returns a pose for RCSB PDB <pdb_code>, also writes this data to
    <pdb_code>.pdb, and optionally calls cleanATOM and/or cleanCRYS

    example:
        pose = pose_from_rcsb("1YY8")
    See also:
        Pose
        pose_from_file
        pose_from_sequence
        load_from_rcsb
        cleanATOM
        cleanCRYS
    """
    pdb_code = pdb_code.upper()
    load_from_rcsb(pdb_code)
    if ATOM:
        cleanATOM(pdb_code + ".pdb")
        pdb_code = pdb_code + ".clean"
    if CRYS:
        cleanCRYS(pdb_code + ".pdb")
        pdb_code = pdb_code + ".mono"
    pose = rosetta.core.import_pose.pose_from_file(pdb_code + ".pdb")
    return pose


# retreives pdbData from rcsb for  <pdb_code>
# ADD NAMING OPTION
def load_fasta_from_rcsb( pdb_code , fasta_outfile ):
    if pdb_code:    # if something input...
        pdb_code = pdb_code.upper()
        try:
            filename = urllib_urlretrieve('http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=' + pdb_code)[0]
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
                if os.path.exists( os.getcwd() + '/' + fasta_outfile ): print( 'the file {} already exists, this file will be overwritten'.format(fasta_outfile) )
                #if input('Do you want to overwrite ' + pdbCode + '.pdb')
                pdb_file = open(fasta_outfile,'w')
                pdb_file.writelines(pdb_data)
                pdb_file.close()

                print( 'PDB {} sequence successfully loaded from rcsb into {}'.format(pdb_code, fasta_outfile) )
            else:
                raise IOError('Invalid PDB code')
        os.remove(filename)    # remove temp file
