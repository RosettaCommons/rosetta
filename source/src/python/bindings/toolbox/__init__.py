# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.
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

# Imports.
import os
import urllib

import rosetta
import rosetta.core.scoring.dssp


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
	p = rosetta.pose_from_pdb(pdbfilename)
	generate_resfile_from_pose(p, resfilename, input_sc)


def mutate_residue(pose, resid, new_res):
	"""
	Replaces the residue at <resid> in <pose> with <new_res>.

	Note: <new_res> is the single letter name for the desired ResidueType.

	Example:
	    mutate_residue(pose, 30, 'A')
	See also:
	    Pose
	    PackRotamersMover
	"""
	if not pose.is_fullatom():
		IOError("mutate_residue only works with full-atom poses.")

	scorefxn = rosetta.create_score_function("talaris2013")
	pack_task = rosetta.TaskFactory.create_packer_task(pose)
	pack_task.initialize_from_command_line()

	v1 = rosetta.utility.vector1_bool()
	mut_res = rosetta.aa_from_oneletter_code(new_res)

	for i in range(1, 21):
		if (i == mut_res):
			v1.append(True)
		else:
			v1.append(False)

	for i in range(1, pose.total_residue() + 1):
		if (i != resid):
			pack_task.nonconst_residue_task(i).prevent_repacking()

	pack_task.nonconst_residue_task(resid).restrict_absent_canonical_aas(v1)

	packer = rosetta.protocols.simple_moves.PackRotamersMover(scorefxn,
	                                                          pack_task)
	packer.apply(pose)
	return pose


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
        pose_from_pdb
        pose_from_rcsb
        cleanCRYS
    """
    if not edit:
        edit = 255
    # why is it pdb_file[:-4] the whole way?
    if os.path.exists(os.getcwd() + '/' + pdb_file):
        print "If the file", pdb_file[:edit] + ".clean.pdb already exists," + \
              "it will be overwritten."
        with open(pdb_file) as f_in:
            with open(pdb_file[:edit] + ".clean.pdb", 'w') as f_out:
                for line in f_in:
                    if line[:4] == "ATOM":
                        f_out.write(line)
        print "PDB", pdb_file, "successfully cleaned, non-ATOM lines removed."
        print "Clean data written to", pdb_file[:edit] + ".clean.pdb."
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
        pose_from_pdb
        pose_from_rcsb
        cleanATOM
    """
    if os.path.exists(os.getcwd() + '/' + pdb_file):
        print "If the file", pdb_file[:-4] + ".mono.pdb already exists, " + \
              "it will be overwritten."
        pose = rosetta.pose_from_pdb(pdb_file)
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
           print "PDB", pdb_file, "successfully cleaned, redundant " + \
                 "monomers removed."
           print "Monomer data written to", pdb_file[:-4] + ".mono.pdb."
        else:
            print pdb_file, "is not a " + str(int(olig)) + "-mer."
    else:
        raise IOError("No such file or directory named " + pdb_file)


def load_from_rcsb(pdb_code, pdb_filename = None):
    """
    Writes PDB data for RCSB data for <pdb_code> into <pdb_filename>. If not
    specified, outputs file to <pdb_code>.pdb.

    Example:
        load_from_rcsb("1YY8")
    See also:
        Pose
        pose_from_pdb
        pose_from_rcsb
        pose_from_sequence
        cleanATOM
        cleanCRYS
    """
    pdb_code = pdb_code.upper()
    try:
        temp = urllib.urlretrieve("http://www.rcsb.org/pdb/files/" +
                                      pdb_code + ".pdb")[0]
    except:
        raise IOError("Cannot access the PDB database, please check your " +
                      " Internet access.")
    else:
        if (os.path.getsize(temp) > 1500):
            # Arbitrarily 1500... else pdb_code was invalid.
            with open(temp) as f:
                pdb_data = f.readlines()

            if not pdb_filename:
                pdb_filename = pdb_code + ".pdb"
            if os.path.exists(os.getcwd() + '/' + pdb_filename):
                print "The file", pdb_filename, "already exists; this file",
                print "will be overwritten."
            with open(pdb_filename, 'w') as f:
                f.writelines(pdb_data)

            print "PDB", pdb_code, "successfully loaded from the RCSB into",
            print pdb_filename + '.'
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
        pose_from_pdb
        pose_from_sequence
        load_from_rcsb
        cleanATOM
        cleanCRYS
    """
    load_from_rcsb(pdb_code)
    if ATOM:
        cleanATOM(pdb_code + ".pdb")
        pdb_code = pdb_code + ".clean"
    if CRYS:
        cleanCRYS(pdb_code + ".pdb")
        pdb_code = pdb_code + ".mono"
    pose = rosetta.pose_from_pdb(pdb_code + ".pdb")
    return pose


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
        pose_from_pdb
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
            print num + '\n' + seq[:page] + '\n' + sec[:page] + '\n'
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
