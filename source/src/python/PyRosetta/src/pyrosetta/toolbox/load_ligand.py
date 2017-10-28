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

from __future__ import print_function

import os
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
    import openbabel
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

        # Add the new ResidueType to the pose
        rts = pose.conformation().modifiable_residue_type_set_for_conf( core.chemical.FULL_ATOM_t )
        rts.add_base_residue_type( name )
        pose.conformation().reset_residue_type_set_for_conf( rts )

        # fill the pose
        pose_from_file( pose , pose.residue_type_set_for_pose() , name + '_0001.pdb')
    else:
    # permanent solution, add to .params list
        add_cid_to_database( cid , name )

        # fill the pose
        pose_from_file( pose , name + '_0001.pdb' )
    return pose

# returns a pose containing a ligand
def pose_from_params( filename , params_list ):
    pose = Pose()

    rts = pose.conformation().modifiable_residue_type_set_for_conf( core.chemical.FULL_ATOM_t )
    rts.read_files_for_base_residue_types( Vector1(params_list) )
    pose.conformation().reset_residue_type_set_for_conf( rts )

    pose_from_file( pose , pose.residue_type_set_for_pose() , filename )
    return pose
