# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import argparse
import os
import pyrosetta
import shutil


database = pyrosetta._rosetta_database_from_env()
if database is None:
    raise RuntimeError("Could not find the PyRosetta database.")
DATABASE_FILES = [
    os.path.join(database, "chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_nonnatural/IGU.params"),
    os.path.join(database, "chemical/residue_type_sets/fa_standard/residue_types/nucleic/dna/GNP.params"),
    os.path.join(database, "chemical/residue_type_sets/fa_standard/residue_types/sidechain_conjugation/CYX.params"),
    os.path.join(database, "chemical/residue_type_sets/fa_standard/patches/nucleic/rna/3prime5prime_methyl_phosphate.txt"),
    # Test files with PDB_ROTAMERS lines containing subfiles
    os.path.join(database, "chemical/residue_type_sets/fa_standard/residue_types/spin_labels/CED.params"),
    os.path.join(database, "chemical/residue_type_sets/fa_standard/residue_types/spin_labels/CED_non_clashing_rotamers.pdb"),
    os.path.join(database, "chemical/residue_type_sets/fa_standard/residue_types/spin_labels/R1A.params"),
    os.path.join(database, "chemical/residue_type_sets/fa_standard/residue_types/spin_labels/R1A_rotamers_nonclashing_moe_saved.pdb"),
    # Test adding a custom PDB_ROTAMERS line containing a subfile
    os.path.join(database, "chemical/residue_type_sets/fa_standard/residue_types/water/TP3.params"),
]

HOH_PDB_ROTAMERS = """
ATOM      1  O   HOH A   1      38.321  30.725  -3.718  1.00 54.44           O  
ATOM      2  H1  HOH A   1      38.569  31.562  -3.313  1.00 54.79           H  
ATOM      3  H2  HOH A   1      38.428  30.954  -4.657  1.00 52.79           H  
TER
ATOM      4  O   HOH A   2      31.363  25.644  21.081  1.00 50.79           O  
ATOM      5  H1  HOH A   2      31.475  26.350  21.755  1.00 50.33           H  
ATOM      6  H2  HOH A   2      30.557  25.199  21.392  1.00 49.91           H  
TER
ATOM      7  O   HOH A   3      38.495  10.813  24.149  1.00 63.20           O  
ATOM      8  H1  HOH A   3      37.729  10.514  23.653  1.00 61.60           H  
ATOM      9  H2  HOH A   3      38.134  11.497  24.721  1.00 62.30           H  
TER
"""


def main(tmp_dir):
    os.chdir(tmp_dir)
    if not pyrosetta.rosetta.basic.was_init_called():
        pyrosetta.init(
            options="-run:constant_seed 1 -ex2",
            extra_options="-out:levels core.init:0 basic.random.init_random_generator:0",
            set_logging_handler="logging",
            silent=True,
        )
    else:
        raise RuntimeError("PyRosetta is already initialized.")

    n_pdb_files = 10

    list_file = os.path.join(tmp_dir, "my_file.list")
    with open(list_file, "w") as f:
        for i in range(n_pdb_files):
            pdb_file = os.path.join(tmp_dir, "listed_{0}.pdb".format(i))
            if i >= n_pdb_files // 2:
                pdb_file += ".gz"
            f.write(pdb_file + os.linesep)
            pose = pyrosetta.pose_from_sequence("A" * (i + 1))
            pyrosetta.dump_pdb(pose, pdb_file)

    for i in range(n_pdb_files):
        pdb_file = os.path.join(tmp_dir, "tmp_{0}.pdb".format(i))
        if i >= n_pdb_files // 2:
            pdb_file += ".gz"
        pose = pyrosetta.pose_from_sequence("G" * (i + 1))
        pyrosetta.dump_pdb(pose, pdb_file)

    for src in DATABASE_FILES:
        if not os.path.isfile(src):
            raise FileNotFoundError("File not found in the PyRosetta database: {0}".format(src))
        dst = os.path.join(tmp_dir, os.path.basename(src))
        shutil.copy(src, dst)
        if dst.endswith("TP3.params"):
            conformers_file = os.path.join(os.path.dirname(dst), f"{os.path.splitext(os.path.basename(dst))[0]}_conformers.pdb")
            with open(dst, "a") as f1, open(conformers_file, "w") as f2:
                f1.write(f"PDB_ROTAMERS {conformers_file}\n")
                f2.write(HOH_PDB_ROTAMERS)


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--tmp_dir', type=str)
    args = parser.parse_args()
    main(args.tmp_dir)
