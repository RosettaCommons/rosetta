# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

"""Brief:   A toolbox of exposed methods for use in PyRosetta.

Remarks: These methods are useful in PyRosetta and are further intended to
         demonstrate proper syntax for various common activities.

         ONLY essential, tested, PyRosetta-exclusive methods should live here!
         Remember, the rest of the Python scripts in "toolbox/" are NOT
         supported, just provided, (ex. environment setup scripts, etc.).

         For those interested, Rosetta has a Surface Area calculator and a
         Radius of Gyration calculator; they are EnergyMethods sa and rg,
         respectively. Create an empty ScoreFunction and use
         ScoreFunction.set_weight to use these calculations.
         Other common geometry calculators are CA_rmsd and all_atom_rmsd.

Author:  Evan Baugh

Edits:   Labonte

"""

from __future__ import print_function

# Imports.
import os

import rosetta
import pyrosetta

from pyrosetta.toolbox.rcsb import load_from_rcsb, pose_from_rcsb, load_fasta_from_rcsb


#Author: Evan Baugh.
#Adding it here, since it has been removed from __init__ for some reason.
# replaces the residue at  <resid>  in  <pose>  with  <new_res>  with repacking
def mutate_residue(pose, mutant_position, mutant_aa,
        pack_radius = 0.0, pack_scorefxn = '' ):
    """
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <pack_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType

    example:
        mutate_residue(pose, 30, A)
    See also:
        Pose
        PackRotamersMover
        MutateResidue
        pose_from_sequence
    """
    #### a MutateResidue Mover exists similar to this except it does not pack
    ####    the area around the mutant residue (no pack_radius feature)
    #mutator = MutateResidue(mutant_position, mutant_aa)
    #mutator.apply(test_pose)

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )


    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = rosetta.core.scoring.get_score_function()

    task = pyrosetta.standard_packer_task(pose)

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = rosetta.utility.vector1_bool()
    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = rosetta.core.chemical.aa_from_oneletter_code(mutant_aa)

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range(1, 21):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == int(mutant_aa) )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task(mutant_position
        ).restrict_absent_canonical_aas(aa_bool)

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    restrict_non_nbrs_from_repacking(pose, mutant_position, task, pack_radius)

    # apply the mutation and pack nearby residues
    #print task
    packer = rosetta.protocols.simple_moves.PackRotamersMover(pack_scorefxn, task)
    packer.apply(pose)

def restrict_non_nbrs_from_repacking(pose, res, task, pack_radius):
    """
    Evan's nbr detection in a function.  Should go in C++
    """

    center = pose.residue( res ).xyz( pose.residue( res ).nbr_atom() )
    print( "Res: pack radius: "+repr(pack_radius) )
    for i in range(1, pose.total_residue() + 1):
        # only pack the mutating residue and any within the pack_radius
            if i == res: continue

            nbr = pose.residue( i ).xyz( pose.residue( i ).nbr_atom() )
            dist = nbr.distance(center)
            if dist > pack_radius:
                task.nonconst_residue_task(i).prevent_repacking()
            else:
                task.nonconst_residue_task(i).restrict_to_repacking()

    #print task
    return task


def generate_resfile_from_pose(pose, resfilename, input_sc = True):
    """
    Writes a resfile for <pose> named <resfilename>, optionally allowing
    input side chains to be used in packing.

    Example:
        generate_resfile_from_pose(pose, "1YY8.resfile")
    See also:
        Pose
        PackRotamersMover
        TaskFactory
    """
    with open(resfilename, 'w') as f:
        id = "NATRO"
        start = ''
        if input_sc:
            start = "USE_INPUT_SC\n"
        f.write(start + "start\n")
        info = pose.pdb_info()
        # pose_from_sequence returns empty PDBInfo; Pose() makes NULL.
        if info and info.nres():
            for i in range(1, pose.total_residue() + 1):
                num = pose.pdb_info().number(i)
                chain = pose.pdb_info().chain(i)
                f.write(str(num).rjust(4) + str(chain).rjust(3) +
                        str(id).rjust(7) + "  \n")
        else:
            for i in range (1, pose.total_residue() + 1):
                num = i
                chain = ' '
                f.write(str(num).rjust(4) + str(chain).rjust(3) +
                        str(id).rjust(7) + "  \n")


def generate_resfile_from_pdb(pdbfilename, resfilename, input_sc = True ):
	"""
	Writes a resfile for PDB file <pdbfilename> named <resfilename>,
	optionally allowing input side chains to be used in packing.

	Example:
	    generate_resfile_from_pdb("1YY8.pdb", "1YY8.resfile")
	See also:
	    Pose
	    PackRotamersMover
	    TaskFactory
	"""
	p = rosetta.core.import_pose.pose_from_file(pdbfilename)
	generate_resfile_from_pose(p, resfilename, input_sc)


def cleanATOM(pdb_file, edit = -4):
    """
    Writes a PDB file from <pdb_file> with all lines not beginning with
    ATOM removed tp <pdb_file>.clean.pdb.

    Note: the second argument, <edit>, is for PDB files not ending in .pdb.

    Example:
        cleanATOM("1YY9.pdb")
    See also:
        Pose
        Pose.dump_pdb
        pose_from_file
        pose_from_rcsb
        cleanCRYS
    """
    if not edit:
        edit = 255
    # why is it pdb_file[:-4] the whole way?
    if os.path.exists(os.getcwd() + '/' + pdb_file):
        print( "If the file", pdb_file[:edit] + ".clean.pdb already exists," + \
              "it will be overwritten." )
        with open(pdb_file) as f_in:
            with open(pdb_file[:edit] + ".clean.pdb", 'w') as f_out:
                for line in f_in:
                    if line[:4] == "ATOM":
                        f_out.write(line)
        print( "PDB", pdb_file, "successfully cleaned, non-ATOM lines removed." )
        print( "Clean data written to", pdb_file[:edit] + ".clean.pdb." )
    else:
        raise IOError("No such file or directory named " + pdb_file)


def cleanCRYS(pdb_file, olig = 2):
    """
    Removes redundant crystal contacts and isolates a monomer by writing a PDB
    file for a monomer of <pdb_file>, if it is an <olig>-mer, to
    <pdb_file>.mono.

    Note: This is by simple sequence comparison.

    Example:
        cleanCRYS("1YY8.pdb", 2)
    See also:
        Pose
        Pose.dump_pdb
        pose_from_file
        pose_from_rcsb
        cleanATOM
    """
    if os.path.exists(os.getcwd() + '/' + pdb_file):
        print( "If the file", pdb_file[:-4] + ".mono.pdb already exists, " + \
              "it will be overwritten." )
        pose = rosetta.core.import_pose.pose_from_file(pdb_file)
        tot = pose.total_residue()
        seq = pose.sequence()
        frags = [''] * olig
        match = [False] * (olig - 1)
        olig = float(olig)
        frac = int(round(tot / olig))
        for f in range(int(olig)):
            frags[f] = seq[:frac]
            seq = seq[frac:]
        for f in range(int(olig-1)):
            match[f] = (frags[0] == frags[f + 1])
        if sum(match) == (olig - 1):
           for i in range(frac * int(olig - 1)):
               pose.delete_polymer_residue(frac + 1)
           pose.dump_pdb(pdb_file[:-4] + ".mono.pdb")
           print( "PDB", pdb_file, "successfully cleaned, redundant " + \
                 "monomers removed." )
           print( "Monomer data written to", pdb_file[:-4] + ".mono.pdb." )
        else:
            print( pdb_file, "is not a " + str(int(olig)) + "-mer." )
    else:
        raise IOError("No such file or directory named " + pdb_file)


def get_secstruct(pose, output=True, space=8, page=80):
    """
    Predicts the secondary structure of <pose>, loading this data into
    the pose's secstruct information and printing the prediction to screen.

    <output> determines if the information is printed to the screen or not;
    <space> and <page> determine formatting.

    Example:
        get_secstruct(pose)
    See also:
        Pose
        Pose.secstruct
        Pose.sequence
        pose_from_file
        pose_from_sequence
        pose_from_rcsb
    """
    dssp = rosetta.core.scoring.dssp.Dssp(pose)
    dssp.insert_ss_into_pose(pose)
    seq = pose.sequence()
    sec = pose.secstruct()
    count = 0
    while len(seq):
        count += 1
        num = ''
        for i in range(page * (count - 1) + 1, page * count, space):
            num += str(i).ljust(space)
        if output:
            print( num + '\n' + seq[:page] + '\n' + sec[:page] + '\n' )
        seq = seq[page:]
        sec = sec[page:]


def get_hbonds(pose):
    """
    Returns an HBondSet of the hydrogen bonding in <pose>.

    Note: More info can be found in rosetta.core.scoring.hbonds.

    Example:
        hbset = get_hbonds(pose)
        hbset.show()
        hbset.hbond(1).acc_res()
    See also:
        Pose
        Energies
        HBondSet
        ScoreFunction
        create_score_function
    """
    pose.update_residue_neighbors()
    hbset = rosetta.core.scoring.hbonds.HBondSet()
    rosetta.core.scoring.hbonds.fill_hbond_set(pose, False, hbset)
    return hbset
