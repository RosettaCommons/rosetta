#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   cleaning.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

import os

from rosetta import Pose
from rosetta import pose_from_file

# removes non ATOM lines from  <pdb_file>  and writes to  <out_file>
def cleanATOM( pdb_file , out_file = '', edit = -4 ):
    """
    Writes all lines in the PDB file  <pdb_file>  beginning with "ATOM" or
    "TER" into  <out_file>  (defaults to  <pdb_file>.clean.pdb)
    note: the third argument, <edit>, if for PDB files not ending in .pdb

    example:
        cleanATOM('1YY9.pdb')
    See also:
        Pose
        Pose.dump_pdb
        pose_from_file
        pose_from_rcsb
        cleanCRYS
    """
    # an optional argument for PDB files not ending in .pdb
    if not edit:
        edit = 255
    # if the file exists
    if os.path.exists( os.getcwd() + '/' + pdb_file ):
        # find all ATOM and TER lines
        fid = open(pdb_file,'r')
        data = fid.readlines()
        fid.close()
        good = []
        for i in data:
            if i[:5] == 'ATOM ' or i[:4] == 'TER ':
                # add your preference rules for ligands, DNA, water, etc.
                good.append(i)
        # default output file to  <pdb_file>.clean.pdb
        if not out_file:
            out_file = pdb_file[:edit]+'.clean.pdb'
        # write the found lines
        print 'if the file',out_file,'already exists, it will be overwritten'
        fid = open(out_file,'w')
        fid.writelines(good)
        fid.close()
        print 'PDB',pdb_file,'successfully cleaned, non-ATOM lines removed\nclean data written to',out_file
        return True
    else:
        print 'No such file or directory named '+pdb_file
        return False

# if you would prefer a simpler call using grep, it looks something like this
#    os.system("grep \"ATOM\" %s.pdb > %s.clean.pdb"%(pdb_file[:edit],pdb_file[:edit]))

# removes redundant crystal contacts, isolate monomer
def cleanCRYS( pdb_file , olig = 2 , out_file = '' ):
    """
    Writes a PDB file for a monomer of  <pdb_file>  if it is a  <olig>-mer
    to  <out_file>  (defaults to  <pdb_file>.mono.pdb)
    note: this is a simple sequence comparison

    example:
        cleanCRYS('1YY8.pdb',2)
    See also:
        Pose
        Pose.dump_pdb
        pose_from_file
        pose_from_rcsb
        cleanATOM
    """
    # if the file exists
    if os.path.exists( os.getcwd() + '/' + pdb_file ):
        # load in the PDB...this is really just to get the sequence
        pose = pose_from_file(pdb_file)
        tot = pose.total_residue()
        seq = pose.sequence()
        # generate sequence fragments until
        frags = ['']*olig
        match = [False]*(olig-1)
        olig = float(olig)
        frac = int(round(tot/olig))
        for f in range(int(olig)):
            frags[f] = seq[:frac]
            seq = seq[frac:]
        # determine if sequence fragments are identical
        for f in range(int(olig-1)):
            match[f] = (frags[0]==frags[f+1])
        # if the protein has repeats, delete all other residues
        if sum(match)==(olig-1):
            for i in range(frac*int(olig-1)):
                pose.delete_polymer_residue(frac+1)    # I hope this works!
            # write the new pdb file
            if not out_file:
                out_file = pdb_file[:-4]+'.mono.pdb'
            print 'if the file',out_file,' already exists, it will be overwritten'
            pose.dump_pdb(out_file)
            print 'PDB',pdb_file,'successfully cleaned, redundant monomers removed\nmonomer data written to',out_file
            return True
        else:
            print pdb_file,'is not a '+str(int(olig))+'-mer'
            return False
    else:
        print 'No such file or directory named '+pdb_file
        return False


