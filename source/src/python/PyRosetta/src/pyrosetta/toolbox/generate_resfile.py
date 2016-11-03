#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   generate_resfile.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

# adapted from original code by Sid Chaudhury

from __future__ import print_function

from pyrosetta import Pose, pose_from_file

# writes a specified resfile for the user, defaults to packing w/ input SC
def generate_resfile_from_pose( pose , resfilename ,
        pack = True , design = False , input_sc = True ,
        freeze = [] , specific = {} ):
    """
    Writes a resfile for  <pose>  named  <resfilename>
       <pack> = True allows packing by default
       <design> = True allows design using all amino acids by default
       <input_sc> = True allows usage of the original side chain conformation
       <freeze> is an optional list of (pose) residue numbers to exclude
            (preserve the side chain conformations of these residues)
       <specific> is an optional dictionary with (pose) residue numbers as keys
            and resfile keywords as corresponding values
            (for setting individual residue options, it may be easier to add
            these numbers to freeze and edit the resfile manually)

    example:
        generate_resfile_from_pose(pose,'1YY8.resfile')
    See also:
        Pose
        PackRotamersMover
        TaskFactory
    """
    # determine the header, default settings
    header = ''
    if pack:
        if not design:
            header += 'NATAA\n'
        else:
            header += 'ALLAA\n# ALLAA will NOT work on bridged Cysteines\n'
    else:
        header += 'NATRO\n'
    if input_sc:
        header += 'USE_INPUT_SC\n'
    to_write = header + 'start\n'
    # add  <freeze>  list to  <specific>  dict
    for i in freeze:
        specific[i] = 'NATRO'
    #  <specific>  is a dictionary with keys() as pose resi numbers
    #    and values as resfile keywords (PIKAA
    # use PDBInfo object to write the resfile
    info = pose.pdb_info()
    # pose_from_sequence returns empty PDBInfo, Pose() makes NULL
    if info and info.nres():
        for i in specific.keys():
            num = pose.pdb_info().number(i)
            chain = pose.pdb_info().chain(i)
            to_write += str(num).rjust(4) + str(chain).rjust(3) + '  ' + specific[i] + '  \n'
    else:
        for i in specific.keys():
            num = i
            chain = ' '
            to_write += str(num).rjust(4) + str(chain).rjust(3) + '  ' + specific[i] + '  \n'
    f = open(resfilename,'w')
    f.write(to_write)
    f.close()

# this is silly, as with cleanCRYS, later implement a Biopy PDBParser method
def generate_resfile_from_pdb( pdbfilename , resfilename ,
        pack = True , design = False , input_sc = True ,
        freeze = [] , specific = {} ):
	"""
    Writes a resfile for the PDB file <pdbfilename>  named  <resfilename>
       <pack> = True allows packing by default
       <design> = True allows design using all amino acids by default
       <input_sc> = True allows usage of the original side chain conformation
       <freeze> is an optional list of (pose) residue numbers to exclude
            (preserve the side chain conformations of these residues)
       <specific> is an optional dictionary with (pose) residue numbers as keys
            and resfile keywords as corresponding values
            (for setting individual residue options, it may be easier to add
            these numbers to freeze and edit the resfile manually)

	example:
	    generate_resfile_from_pdb('1YY8.pdb','1YY8.resfile')
	See also:
	    generate_resfile_from_pose
	    Pose
	    PackRotamersMover
	    TaskFactory
	"""
	p = pose_from_file(pdbfilename)
	generate_resfile_from_pose(p,resfilename,pack,design,input_sc,freeze,specific)

# this class currently supports NOTHING more than the options above...but that will change! soon?
class ResfileWriter():
    def __init__( self , pose , resfilename , pack = True , design = False , input_sc = True , freeze = [] , specific = {} ):
        self.pose = pose
        self.resfilename = resfilename
        self.pack = pack
        self.design = design
        self.input_sc = input_sc
        self.freeze = freeze
        self.specific = specific

    def write_resfile( self , resfilename = '' ):
        if not resfilename:
            resfilename = self.resfilename
        generate_resfile_from_pose( self.pose , resfilename , self.pack , self.design , self.input_sc , self.freeze , self.specific )

# make it actually perform from PDB, use Bio.PDBParser
