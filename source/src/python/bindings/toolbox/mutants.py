#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   mutants.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

# mutate_residue is adapted from an original script by Sid Chaudhury

# UNFINISHED...below!

import random
from rosetta import Pose
from rosetta import make_pose_from_sequence
from rosetta import create_score_function
from rosetta import TaskFactory
from rosetta.utility import vector1_bool
from rosetta import aa_from_oneletter_code
from rosetta import PackRotamersMover
from rosetta.core.pose import PDBInfo

# a different version of mutate_residue is provided in PyRosetta v2.0 and
#    earlier that does not optionally repack nearby residues

# replaces the residue at  <resid>  in  <pose>  with  <new_res>  with repacking
def mutate_residue( pose , mutant_position , mutant_aa ,
        pack_radius = 0.0 , pack_scorefxn = '' ):
    """
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <pack_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType

    example:
        mutate_residue(pose,30,A)
    See also:
        Pose
        PackRotamersMover
        MutateResidue
        pose_from_sequence
    """
    #### a MutateResidue Mover exists similar to this except it does not pack
    ####    the area around the mutant residue (no pack_radius feature)
    #mutator = MutateResidue( mutant_position , mutant_aa )
    #mutator.apply( test_pose )

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = Pose()
    test_pose.assign( pose )

    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = create_score_function( 'standard' )

    task = standard_packer_task( test_pose )

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = rosetta.utility.vector1_bool()
    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = aa_from_oneletter_code( mutant_aa )

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range( 1 , 21 ):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == mutant_aa )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task( mutant_position
        ).restrict_absent_canonical_aas( aa_bool )

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue( mutant_position ).nbr_atom_xyz()
    for i in range( 1 , pose.total_residue() + 1 ):
        # only pack the mutating residue and any within the pack_radius
        if not i == mutant_position or center.distance_squared(
                test_pose.residue( i ).nbr_atom_xyz() ) > pack_radius**2:
            task.nonconst_residue_task( i ).prevent_repacking()

    # apply the mutation and pack nearby residues
    packer = PackRotamersMover( pack_scorefxn , task )
    packer.apply( test_pose )

    return test_pose

# returns a pose made from seq, all phi/psi/omega set to 180
def pose_from_sequence( seq , res_type = 'fa_standard' , name = '' , chain_id = 'A' ):
    """
    Returns a pose generated from amino acid single letters in  <seq>  using
    the  <res_type>  ResidueType, the new pose's PDBInfo is named  <name>
    and all residues have chain ID  <chain_id>

    example:
        pose=pose_from_sequence('LIGAND')
    See also:
        Pose
        make_pose_from_sequence
        pose_from_file
        pose_from_rcsb
    """
    pose=Pose()
    make_pose_from_sequence(pose,seq,res_type)
    #pdb_info = rosetta.core.pose.PDBInfo(pose.total_residue())    # actual, for other code
    pdb_info = PDBInfo(pose.total_residue())    # create a PDBInfo object
    for i in range(0,pose.total_residue()):
        if pose.residue(i+1).is_protein():
            # set to a more reasonable default
            pose.set_phi(i+1,-150)
            pose.set_psi(i+1,150)
            pose.set_omega(i+1,180)
            # set PDBInfo info for chain and number
        #pdb_info.chain(i+1,chain_id)
        #pdb_info.number(i+1,i+1)
    #### you can alternatively use the deprecated method set_extended_torsions
    ####    which requires a Pose and a Loop object...so make a large loop
    #set_extended_torsions( pose , Loop ( 1 , pose.total_residue() ) )
    # set the PDBInfo
    pose.pdb_info(pdb_info)
    # default name to first 3 letters
    if not name:
        name = seq[:4]
    pose.pdb_info().name(name)
#	print pose
    return pose

# for use with random_sequence
protein_letters = 'ACDEFGHIKLMNPQRSTVWY'
nucleic_letters = ['A[ADE]','G[GUA]','C[CYT]','T[THY]']

# returns a sequence of random protein letters
def random_sequence( length = 1 , letters = protein_letters ):
    return ''.join( [letters[random.randint(0,len(letters)-1)] for i in range( length ) ] )

# the following method may not be useulf to many...but whatever

# display sequence differences
def compare_mutants__( pose1 , pose2 ):
    seq1 = pose1.sequence()
    seq2 = pose2.sequence()
    for i in range(pose1.total_residue()):
        if not seq1[i]==seq2[i]:
            print 'mutation '+seq1[i]+str(i+1)+seq2[i]

# use to compare sequences and sec_struct
def compare_sequences( seq1 , seq2 ):
    for i in range(len(seq1)):
        if not seq1[i]==seq2[i]:
            print 'discrepency at '+seq1[i]+str(i+1)+seq2[i]

# UNFINISHED BELOW HERE!

# look for unmatched hbonds
def compare_hbonds( pose1 , pose2 , Ethresh = .5 , display = False ):
    hblist1 = hbond_list(pose1)
    hbset1 = get_hbonds(pose1)
    hblist2 = hbond_list(pose2)
    hbset2 = get_hbonds(pose2)
    if display:
        pymol = PyMOL_Mover()
    delEn = []
    for i in range(len(hblist1)):
        for j in range(len(hblist2)):
            if hblist1[i][1]==hblist2[j][1] and hblist1[i][2]==hblist2[j][2]:
                # if the energy changes by ~25% or more...
                if abs((hblist2[j][3]-hblist1[i][3])/hblist1[i][3])>Ethresh:
                    # compare bb sc ?
                    delEn.append(hblist1[i])
                    delEn[len(delEn)-1].append(hblist2[j][3])
                    delEn[len(delEn)-1].append(hblist2[j][3]>hblist1[i][3])
                    if display:
                        delEn[len(delEn)-1].append(hbset1.hbond(i))
                hblist1[i] = (0,0,0)
                hblist2[j] = (0,0,0)
                break

    # remove all the matches now
    for i in range(hblist1.count((0,0,0))):
        hblist1.remove((0,0,0))
    for j in range(hblist2.count((0,0,0))):
        hblist2.remove((0,0,0))

    # print the explicit data
    names = [pose1.pdb_info().name(),pose2.pdb_info().name()]
   # is display:    #set this up to avoid lots of ifs

    if len(hblist2)>0:
        print 'hbonds not in',names[0]
    for j in range(len(hblist2)):
        # these hbonds are in pose2 but NOT in pose1
        print hblist2[j][1],'\t',hblist2[j][2]
        if display:
            # send hbonds gained, green
            # align them first!!!!
            hbond = hbset2.hbond(hblist2[j][0])
            donor = pose2.residue(hbond.don_res()).xyz(hbond.don_hatm())
            acptr = pose2.residue(hbond.acc_res()).xyz(hbond.acc_atm())
            pymol.send_point(names[0]+'_gains'+str(j),'green',donor[0],donor[1],donor[2],False,False,'')
            pymol.send_point(names[0]+'_gains'+str(j),'green',acptr[0],acptr[1],acptr[2],False,False,'')

    if len(hblist1)>0:
        print 'hbonds not in',names[1]
    for i in range(len(hblist1)):
        # these hbonds are in pose1 but NOT in pose2
        print hblist1[i][1],'\t',hblist1[i][2]
        if display:
            # send hbonds lost, red
            hbond = hbset1.hbond(hblist1[i][0])
            donor = pose1.residue(hbond.don_res()).xyz(hbond.don_hatm())
            acptr = pose1.residue(hbond.acc_res()).xyz(hbond.acc_atm())
            pymol.send_point(names[0]+'_loses'+str(i),'red',donor[0],donor[1],donor[2],False,False,'')
            pymol.send_point(names[0]+'_loses'+str(i),'red',acptr[0],acptr[1],acptr[2],False,False,'')

    if len(delEn)>0:
        print 'Energy change from',names[0],'to',names[1]
    for k in range(len(delEn)):
        sign = '-'
        color = 'blue'
        group = '_decreases'
        if delEn[k][5]:
            sign = '+'
            color = 'yellow'
            group = '_increases'
        print sign,'in',str(delEn[k][1]).ljust(20),str(delEn[k][2]).ljust(20),'from %.4f to %.4f'%(delEn[k][3],delEn[k][4])
        if display:
            # send energy change, if sign=='+' clr yellow, if sign=='-' clr blue
            hbond = delEn[k][6]
            donor = pose1.residue(hbond.don_res()).xyz(hbond.don_hatm())
            acptr = pose1.residue(hbond.acc_res()).xyz(hbond.acc_atm())
            pymol.send_point(names[0]+group+str(k),color,donor[0],donor[1],donor[2],False,False,'')
            pymol.send_point(names[0]+group+str(k),color,acptr[0],acptr[1],acptr[2],False,False,'')


