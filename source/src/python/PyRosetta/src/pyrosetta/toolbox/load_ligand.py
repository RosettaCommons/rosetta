#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   load_ligand.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

# generate_nonstandard_residue_set is adapted from an original script by Sid Chaudhury

from __future__ import print_function

import os
# its in python!
import openbabel
import molfile_to_params
from pyrosetta import Pose, pose_from_file, init

################################################################################
# methods for obtaining ligand chemical files and producing params .files

# 1. get the ligand file, as .sdf
# retreives sdfdata from rcsb...currently inefficient
def load_from_pubchem( cid , sdffilename = '' ):
    cid = str(cid)
    try:
        filename = urllib.urlretrieve('http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=' + cid + '&disopt=3DDisplaySDF')[0]
    except:
        raise IOError('Cannot access pubchem database, please check your Internet access')
    else:
        if os.path.getsize(filename)>0:    # arbitrary...then cid was invalid
            sdfile = open(filename)
            sdfdata = sdfile.read()
            sdfile.close()

            # default naming
            if not sdffilename:
                sdffilename = 'CID_' + cid + '.sdf'
            sdfile = open(sdffilename,'w')
            sdfile.write(sdfdata)
            sdfile.close()

            print( 'CID {} successfully written to .sdf file {}'.fortma(cid, sdffilename) )

            #mlf = mdl_molfile.read_mdl_molfile(cid+'.sdf')
        else:
            raise IOError('Invalid CID code')
        os.remove(filename)    # remove tmp file

# 2. generate conformers! UNSUPPORTED
#def conformers_from_sdf( sdffilename ):
#    # returns a list of conformers generated from the input .sdf file
#    print 'No conformer generation supported, sry we fail HARD at this'
#    return [sdffilename]

# 3, convert to .mol or ,mdl
# uses openbabel to convert sdf to mdl
def sdf2mdl( sdfile , mdlfilename ):
    if os.path.exists( os.getcwd() + '/' + sdfile ):
        converter = openbabel.OBConversion()
        converter.SetInAndOutFormats('sdf','mdl')
        mol = openbabel.OBMol()
        converter.ReadFile(mol,sdfile)
        print( 'if the file {} already exists, it will be overwritten...'.format(mdlfilename) )
#        os.system("babel  %s %s.pdb"%(sdfile,sdfile[:-4]))
        converted = converter.WriteFile(mol,mdlfilename)
        if converted:
            print( '.mdl file {} successfully written'.format(mdlfilename) )
        else:
            print( 'Conversion Failed! could not produce the .mdl file {} '.format(mdlfilename) )
    else:
        raise IOError('No such file or directory named '+sdfile)

# 4. convert .mdl to .params
# quick wrapper for molfile_to_params with basic functionality
def molfile2params_quick( mdlfile , name ):
    molfile_to_params.main([mdlfile,'-n'+name])

################################################################################
# Temporary solution, load the ligand for this session

# a method for producing non-standard ResidueTypeSets
#    there is a custom (PyRosetta only) method in PyRosetta v2.0beta named
#    generate_nonstandard_residue_set which will work there, but not for
#    newer versions, this method supports older and newer versions of PyRosetta
def generate_nonstandard_residue_set( params_list ):
    """
    Returns a "custom" ResidueTypeSet with the normal ResidueTypes and any
        new ones added as a Vector1 of .params filenames,
        the input  <params_list>

    example(s):
        res_set = generate_nonstandard_residue_set( Vector1( ['ATP.params'] ) )
    See Also:
        Pose
        Residue
        ResidueType
        ResidueTypeSet
    """
    res_set = ChemicalManager.get_instance().nonconst_residue_type_set(
        'fa_standard' )
    atoms = ChemicalManager.get_instance().atom_type_set( 'fa_standard' )
    mm_atoms = ChemicalManager.get_instance().mm_atom_type_set( 'fa_standard' )
    orbitals = ChemicalManager.get_instance().orbital_type_set( 'fa_standard' )
    try:
        # then this PyRosetta is a newer version, sorry, element_sets were added
        #    to the chemical database and changed the syntax of read_files
        elements = ChemicalManager.get_instance().element_set( 'default' )
        res_set.read_files( params_list , atoms , elements , mm_atoms , orbitals )
    except:
        # then this PyRosetta is v2.0 beta or earlier, as this is being written,
        #    we support v2.0 beta, notice the subtle difference below
        res_set.read_files( params_list , atoms , mm_atoms , orbitals )
    return res_set

# perform the above steps in one function call
def params_from_pubchem( cid , name ):
    filename = 'CID_' + str(cid)

    # load from pubchem
    load_from_pubchem( cid )

    # generate conformers
    #files = conformers_from_sdf( filename + '.sdf' )

    # convert to .mdl
    sdf2mdl( filename + '.sdf' , filename + '.mdl' )

    # produce .params
    molfile2params_quick( filename + '.mdl' . name )


################################################################################
# Permanent solution, add the .params to the minirosetta_database of PyRosetta

# this may not work if you manipulated your path variables
database = os.path.abspath( os.environ['PYROSETTA_DATABASE'] )
fa_standard = database + '/chemical/residue_type_sets/fa_standard/'
fa_custom = 'residue_types/custom'

def add_cid_to_database( cid , name ):
    # change directory to the custom database
    start_dir = os.getcwd()
    # if it does not exist
    if not os.path( database + fa_standard + fa_custom ):
        # make a "custom" directory
        os.chdir( database + fa_standard + 'residue_types' )
        os.mkdir( 'custom' )
        # edit residue_type_sets.txt
        os.chdir( database + fa_standard )
        f = open( 'residue_type_sets.txt' , 'w' )
        data = f.readlines()
        data.append( '\n## Custom\n' )
        f.readlines( data )
        f.close()
    os.chdir( database + fa_standard + fa_custom )

    # get the ligand
    params_from_pubchem( cid , name )

    # add the ligand to residue_type_sets.txt
    os.chdir( database + fa_standard )
    f = open( 'residue_type_sets.txt' , 'w' )
    data = f.readlines()
    data.append( fa_custom + '/' + name + '.params\n' )
    f.readlines( data )
    f.close()

    # return to original dir
    os.chdir( start )

    # reinitialize
    init()

################################################################################
# returns a pose of the molecule
def pose_from_pubchem( cid , name , temporary = True ):
    pose = Pose()
    if temporary:
    # the temporary solution, create an ephemeral ResidueSet
        params_from_pubchem( cid , name )

        # generate ResidueSet
        res_set = generate_nonstandard_residue_set( [name] )

        # fill the pose
        pose_from_file( pose , res_set , name + '_0001.pdb')
    else:
    # permanent solution, add to .params list
        add_cid_to_database( cid , name )

        # fill the pose
        pose_from_file( pose , name + '_0001.pdb' )
    return pose

# returns a pose containing a ligand
def pose_from_params( filename , params_list ):
    res_set = generate_nonstandard_residue_set( params_list )
    pose = Pose()
    pose_from_file( pose , res_set , filename )
    return pose
