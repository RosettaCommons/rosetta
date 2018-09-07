import pyrosetta.rosetta as rosetta

# By Michael Pacella
def etable_atom_pair_energies(res1, atom_index_1, res2, atom_index_2, sfxn):
    """Compute the energy of two atoms and return the LJ, solvation and electrostatic
    terms.

    Args:
        res1 (pyrosetta.rosetta.core.conformation.Residue): the residue that contains the
            first atom of interest.
        atom_index_1 (int): index of the desired atom in residue 1
        res2 (pyrosetta.rosetta.core.conformation.Residue): the residue that contains the
            second atom of interest.
        atom_index_2 (int): index of the desired atom in residue 2

    Returns:
        tuple: values of the lj_atr, lj_rep, fa_solv, and fa_elec potentials.

    Usage: lj_atr, lj_rep, solv=etable_atom_pair_energies(res1, atom_index_1, res2, atom_index_2, sfxn)
        Description: given a pair of atoms (specified using a pair of residue objects and
        atom indices) and scorefunction, use the precomputed 'etable' to return
        LJ attractive, LJ repulsive, and LK solvation energies
    """
    score_manager = rosetta.core.scoring.ScoringManager.get_instance()
    etable_ptr = score_manager.etable(sfxn.energy_method_options().etable_type())
    etable = etable_ptr.lock()
    etable_energy = rosetta.core.scoring.etable.AnalyticEtableEnergy(
        etable, sfxn.energy_method_options()
    )

    # Construct coulomb class for calculating fa_elec energies
    coulomb = rosetta.core.scoring.etable.coulomb.Coulomb(sfxn.energy_method_options())

    # Construct AtomPairEnergy container to hold computed energies.
    ape = rosetta.core.scoring.etable.AtomPairEnergy()

    # Set all energies in the AtomPairEnergy to zero prior to calculation.
    ape.attractive, ape.bead_bead_interaction, ape.repulsive, ape.solvation = [0.] * 4

    # Calculate the distance squared and set it in the AtomPairEnergy.
    ape.distance_squared = res1.xyz(atom_index_1).distance_squared(
        res2.xyz(atom_index_2)
    )

    # Evaluate energies from pre-calculated etable, using a weight of 1.0
    # in order to match the raw energies from eval_ci_2b.
    atom1 = res1.atom(atom_index_1)
    atom2 = res2.atom(atom_index_2)
    etable_energy.atom_pair_energy(atom1, atom2, 1.0, ape)

    # Calculate atom-atom scores.
    lj_atr = ape.attractive
    lj_rep = ape.repulsive
    solv = ape.solvation
    fa_elec = coulomb.eval_atom_atom_fa_elecE(
        res1.xyz(atom_index_1),
        res1.atomic_charge(atom_index_1),
        res2.xyz(atom_index_2),
        res2.atomic_charge(atom_index_2),
    )

    return lj_atr, lj_rep, solv, fa_elec


def _atom_pair_energy_table(sfxn, score_type, residue_1, residue_2, threshold=0.):
    """Compute the pairwise atom pair energies for two residues and yield the result
    for each atom on residue 1.

    Args:
        sfxn (pyrosetta.ScoreFunction): the scorefunction to use
        score_type (str): the term of the scorefunction to yield
        residue_1 (pyrosetta.rosetta.core.conformation.Residue): residue 1
        residue_2 (pyrosetta.rosetta.core.conformation.Residue): residue 2
        threshold (float): Minimum absolute value of a score term to include.
            Defaults to 0.

    Yields:
        list: a list of tuples containing the atom index of residue 2 and the value
        of the score term representing all interatomic interactions of a particular
        atom on residue 1.
    """
    terms = ["fa_atr", "fa_rep", "fa_sol", "fa_elec"]

    if score_type not in terms:
        print(
            "please enter a valid score_type: {}, or {}".format(
                ", ".join(terms[:-1]), terms[-1]
            )
        )
        return

    for i in range(1, residue_1.natoms() + 1):
        list_of_interactions = [residue_1.atom_name(i)]
        for j in range(1, residue_2.natoms() + 1):
            score_dict = dict(
                zip(terms, etable_atom_pair_energies(residue_1, i, residue_2, j, sfxn))
            )
            # score_type is guaranteed to be a valid key because we already checked
            list_of_interactions.append(((i, j), score_dict[score_type]))
        yield list_of_interactions


def dump_atom_pair_energy_table(
    sfxn, score_type, residue_1, residue_2, output_filename
):
    """Compute of all pairwise atom pair energies for the complete list of atoms
    contained by residue_1 and residue_2 using a specified score_type in the provided
    sfxn and write the results to a csv file.

    Notes:
        Only non-zero scores are written

    Args:
        sfxn (pyrosetta.ScoreFunction): the scorefunction to use
        score_type (str): the term of the scorefunction to consider
        residue_1 (pyrosetta.rosetta.core.conformation.Residue): residue 1
        residue_2 (pyrosetta.rosetta.core.conformation.Residue): residue 2
        output_filename (str): name of the file to which the results will be written
    """

    list_of_res1_atoms = list(
        _atom_pair_energy_table(sfxn, score_type, residue_1, residue_2)
    )
    list_of_res1_atoms = [l[-1] for l in list_of_res1_atoms]

    with open(output_filename + "_" + str(score_type) + ".csv", "wb") as f:
        writer = csv.writer(f)
        writer.writerows(list_of_res1_atoms)


def print_atom_pair_energy_table(sfxn, score_type, residue_1, residue_2, threshold):
    """Compute of all pairwise atom pair energies for the complete list of atoms
    contained by residue_1 and residue_2 using a specified score_type in the provided
    sfxn and print the values that exceed some threshold value.

    Args:
        sfxn (pyrosetta.ScoreFunction): the scorefunction to use
        score_type (str): the term of the scorefunction to consider
        residue_1 (pyrosetta.rosetta.core.conformation.Residue): residue 1
        residue_2 (pyrosetta.rosetta.core.conformation.Residue): residue 2
        threshold (float): Minimum absolute value of a score term to include.
            Defaults to 0.
    """
    for (res1_atm_idx, res2_atm_idx), score in _atom_pair_energy_table(
        sfxn, score_type, residue_1, residue_2, threshold
    ):
        print(
            "{0} {1} {2}".format(
                residue_1.atom_name(res1_atm_idx),
                residue_2.atom_name(res2_atm_idx),
                score,
            )
        )


def _reisude_pair_energies(res, pose, sfxn, score_type, threshold=0.):
    """Compute the scores of a particular residue with all other residues in a Pose
    and yield them.

    Args:
        res (int): Pose-numbered index of the residue of interest
        pose (pyrosetta.Pose): the target Pose
        sfxn (pyrosetta.ScoreFunction): the scorefunction to use
        score_type (str): the term of the scorefunction to return
        threshold (float): Minimum absolute value of a score term to include.
            Defaults to 0.

    Yields:
        tuple: (i, score) where i is a Pose-numbered residue and score is the score of
            residue i with the user-specified residue, res.
    """
    sfxn.score(pose)
    for i in range(1, pose.total_residue() + 1):
        emap = rosetta.core.scoring.EMapVector()
        sfxn.eval_ci_2b(pose.residue(res), pose.residue(i), pose, emap)
        if abs(emap[score_type]) > threshold:
            yield (i, emap[score_type])


def dump_residue_pair_energies(res, pose, sfxn, score_type, output_filename):
    """Compute the scores of a particular residue with all other residues in a Pose
    and write the results to a csv file.

    Notes:
        Only non-zero scores are returned

    Args:
        res (int): Pose-numbered index of the residue of interest
        pose (pyrosetta.Pose): the target Pose
        sfxn (pyrosetta.ScoreFunction): the scorefunction to use
        score_type (str): the term of the scorefunction to return
        output_filename (str): name of the file to which the results will be written
    """
    list_of_scores = list(_reisude_pair_energies(res, pose, sfxn, score_type))

    with open(output_filename + "_" + str(score_type) + ".csv", "wb") as f:
        writer = csv.writer(f)
        writer.writerows(list_of_scores)


def print_residue_pair_energies(res, pose, sfxn, score_type, threshold):
    """Compute the scores of a particular residue with all other residues in a Pose
    and print the values that exceed some threshold value.

    Args:
        res (int): Pose-numbered index of the residue of interest
        pose (pyrosetta.Pose): the target Pose
        sfxn (pyrosetta.ScoreFunction): the scorefunction to use
        score_type (str): the term of the scorefunction to return
        threshold (float): Minimum absolute value of a score term to include.
            Defaults to 0.
    """
    for res_idx, score in _reisude_pair_energies(
        res, pose, sfxn, score_type, threshold
    ):
        print("{0} {1}  {2}".format(pose.residue(res_idx).name1(), res_idx, score))
