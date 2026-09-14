#!/usr/bin/env python
# :noTabs=true:
# -*- coding: utf-8 -*-

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   BindCraft.py
## @brief  Tests for the PyRosetta function of BindCraft
## @brief       (https://github.com/martinpacesa/BindCraft),
## @brief  mainly to make sure we don't inadvertently break them.

'''What follows is the contents of https://github.com/martinpacesa/BindCraft/blob/main/functions/pyrosetta_utils.py
Reproduced here under the MIT license:

    MIT License

    Copyright (c) 2024 Martin Pacesa

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

Portions which aren't strictly testing the PyRosetta functionality have been commented out with `##%##`
'''
###################################################################################################################
####################################
################ PyRosetta functions
####################################
### Import dependencies
import os
import pyrosetta as pr
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.io import pose_from_pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
##%##from .generic_utils import clean_pdb
##%##from .biopython_utils import hotspot_residues

# Rosetta interface scores
def score_interface(pdb_file, binder_chain="B"):
    # load pose
    pose = pr.pose_from_pdb(pdb_file)

    # analyze interface statistics
    iam = InterfaceAnalyzerMover()
    interface = "A_B"
    docking_partners_type = getattr(pr.rosetta.core.pose, "DockingPartners", None)
    if docking_partners_type is not None:
        interface = docking_partners_type.docking_partners_from_string(interface)
    iam.set_interface(interface)
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)

##%## # All this is commented out to avoid hotspot_residues() usage
##%##    # Initialize dictionary with all amino acids
##%##    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
##%##
##%##    # Initialize list to store PDB residue IDs at the interface
##%##    interface_residues_set = hotspot_residues(pdb_file, binder_chain)
##%##    interface_residues_pdb_ids = []
##%##
##%##    # Iterate over the interface residues
##%##    for pdb_res_num, aa_type in interface_residues_set.items():
##%##        # Increase the count for this amino acid type
##%##        interface_AA[aa_type] += 1
##%##
##%##        # Append the binder_chain and the PDB residue number to the list
##%##        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")
##%##
##%##    # count interface residues
##%##    interface_nres = len(interface_residues_pdb_ids)
    interface_nres = 50 ##%## Semi-arbitrary to just skip the interface issue

##%##    # Convert the list into a comma-separated string
##%##    interface_residues_pdb_ids_str = ','.join(interface_residues_pdb_ids)
##%##
##%##    # Calculate the percentage of hydrophobic residues at the interface of the binder
##%##    hydrophobic_aa = set('ACFILMPVWY')
##%##    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
##%##    if interface_nres != 0:
##%##        interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100
##%##    else:
##%##        interface_hydrophobicity = 0

    # retrieve statistics
    interfacescore = iam.get_all_data()
    interface_sc = interfacescore.sc_value # shape complementarity
    interface_interface_hbonds = interfacescore.interface_hbonds # number of interface H-bonds
    interface_dG = iam.get_interface_dG() # interface dG
    interface_dSASA = iam.get_interface_delta_sasa() # interface dSASA (interface surface area)
    interface_packstat = iam.get_interface_packstat() # interface pack stat score
    interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100 # ratio of dG/dSASA (normalised energy for interface area size)
    buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="0" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
##%## # Unfortunately, we don't have dalphaball on the test server, so we turn it off.
##%##    buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
    interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    if interface_nres != 0:
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100 # Hbonds per interface size percentage
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100 # Unsaturated H-bonds per percentage
    else:
        interface_hbond_percentage = None
        interface_bunsch_percentage = None

    # calculate binder energy score
    chain_design = ChainSelector(binder_chain)
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_design)
    binder_score = tem.calculate(pose)

    # calculate binder SASA fraction
    bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(chain_design)
    binder_sasa = bsasa.calculate(pose)

    if binder_sasa > 0:
        interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
    else:
        interface_binder_fraction = 0

    # calculate surface hydrophobicity
    binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}[binder_chain]

    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core = False, pick_boundary = False, pick_surface = True)
    surface_res = layer_sel.apply(binder_pose)

    exp_apol_count = 0
    total_count = 0

    # count apolar and aromatic residues at the surface
    for i in range(1, len(surface_res) + 1):
        if surface_res[i] == True:
            res = binder_pose.residue(i)

            # count apolar and aromatic residues as hydrophobic
            if res.is_apolar() == True or res.name() == 'PHE' or res.name() == 'TRP' or res.name() == 'TYR':
                exp_apol_count += 1
            total_count += 1

    surface_hydrophobicity = exp_apol_count/total_count

    # output interface score array and amino acid counts at the interface
    interface_scores = {
    'binder_score': binder_score,
    'surface_hydrophobicity': surface_hydrophobicity,
    'interface_sc': interface_sc,
    'interface_packstat': interface_packstat,
    'interface_dG': interface_dG,
    'interface_dSASA': interface_dSASA,
    'interface_dG_SASA_ratio': interface_dG_SASA_ratio,
    'interface_fraction': interface_binder_fraction,
##%##    'interface_hydrophobicity': interface_hydrophobicity,
    'interface_nres': interface_nres,
    'interface_interface_hbonds': interface_interface_hbonds,
    'interface_hbond_percentage': interface_hbond_percentage,
    'interface_delta_unsat_hbonds': interface_delta_unsat_hbonds,
    'interface_delta_unsat_hbonds_percentage': interface_bunsch_percentage
    }

    # round to two decimal places
    interface_scores = {k: round(v, 2) if isinstance(v, float) else v for k, v in interface_scores.items()}

    return interface_scores
##%##    return interface_scores, interface_AA, interface_residues_pdb_ids_str

# align pdbs to have same orientation
def align_pdbs(reference_pdb, align_pdb, reference_chain_id, align_chain_id):
    # initiate poses
    reference_pose = pr.pose_from_pdb(reference_pdb)
    align_pose = pr.pose_from_pdb(align_pdb)

    align = AlignChainMover()
    align.pose(reference_pose)

    # If the chain IDs contain commas, split them and only take the first value
    reference_chain_id = reference_chain_id.split(',')[0]
    align_chain_id = align_chain_id.split(',')[0]

    # Get the chain number corresponding to the chain ID in the poses
    reference_chain = pr.rosetta.core.pose.get_chain_id_from_chain(reference_chain_id, reference_pose)
    align_chain = pr.rosetta.core.pose.get_chain_id_from_chain(align_chain_id, align_pose)

    align.source_chain(align_chain)
    align.target_chain(reference_chain)
    align.apply(align_pose)

    # Overwrite aligned pdb
    align_pose.dump_pdb(align_pdb)
##%##    clean_pdb(align_pdb)

# calculate the rmsd without alignment
def unaligned_rmsd(reference_pdb, align_pdb, reference_chain_id, align_chain_id):
    reference_pose = pr.pose_from_pdb(reference_pdb)
    align_pose = pr.pose_from_pdb(align_pdb)

    # Define chain selectors for the reference and align chains
    reference_chain_selector = ChainSelector(reference_chain_id)
    align_chain_selector = ChainSelector(align_chain_id)

    # Apply selectors to get residue subsets
    reference_chain_subset = reference_chain_selector.apply(reference_pose)
    align_chain_subset = align_chain_selector.apply(align_pose)

    # Convert subsets to residue index vectors
    reference_residue_indices = get_residues_from_subset(reference_chain_subset)
    align_residue_indices = get_residues_from_subset(align_chain_subset)

    # Create empty subposes
    reference_chain_pose = pr.Pose()
    align_chain_pose = pr.Pose()

    # Fill subposes
    pose_from_pose(reference_chain_pose, reference_pose, reference_residue_indices)
    pose_from_pose(align_chain_pose, align_pose, align_residue_indices)

    # Calculate RMSD using the RMSDMetric
    rmsd_metric = RMSDMetric()
    rmsd_metric.set_comparison_pose(reference_chain_pose)
    rmsd = rmsd_metric.calculate(align_chain_pose)

    return round(rmsd, 2)

# Relax designed structure
def pr_relax(pdb_file, relaxed_pdb_path):
    if not os.path.exists(relaxed_pdb_path):
        # Generate pose
        pose = pr.pose_from_pdb(pdb_file)
        start_pose = pose.clone()

        ### Generate movemaps
        mmf = MoveMap()
        mmf.set_chi(True) # enable sidechain movement
        mmf.set_bb(True) # enable backbone movement, can be disabled to increase speed by 30% but makes metrics look worse on average
        mmf.set_jump(False) # disable whole chain movement

        # Run FastRelax
        fastrelax = FastRelax()
        scorefxn = pr.get_fa_scorefxn()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf) # set MoveMap
        fastrelax.max_iter(200) # default iterations is 2500
        fastrelax.min_type("lbfgs_armijo_nonmonotone")
        fastrelax.constrain_relax_to_start_coords(True)
        fastrelax.apply(pose)

        # Align relaxed structure to original trajectory
        align = AlignChainMover()
        align.source_chain(0)
        align.target_chain(0)
        align.pose(start_pose)
        align.apply(pose)

        # Copy B factors from start_pose to pose
        for resid in range(1, pose.total_residue() + 1):
            if pose.residue(resid).is_protein():
                # Get the B factor of the first heavy atom in the residue
                bfactor = start_pose.pdb_info().bfactor(resid, 1)
                for atom_id in range(1, pose.residue(resid).natoms() + 1):
                    pose.pdb_info().bfactor(resid, atom_id, bfactor)

        # output relaxed and aligned PDB
        pose.dump_pdb(relaxed_pdb_path)
##%##        clean_pdb(relaxed_pdb_path)

###################################################################################################################
# End transclusion

import os
import unittest
import shutil
import pyrosetta
from pyrosetta import rosetta
import tempfile

class BindCraftTest(unittest.TestCase):
    '''These aren't intended as extensive tests of functionality, just white-box smoke tests for proper interface functioning.'''
    @classmethod
    def setUpClass(cls):
        # These options match those from the init() call of ./bindcraft.py (modulo -constant_seed, -mute and -holes:dalphaball)
        pyrosetta.init(extra_options = "-constant_seed -ignore_unrecognized_res -ignore_zero_occupancy -corrections::beta_nov16 true -relax:default_repeats 1")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
        print( pyrosetta.version() )
        cls.workdir = tempfile.TemporaryDirectory()
        os.makedirs('.test.output', exist_ok=True) # In case it doesn't exist for local testings, etc.
        os.chdir('.test.output')

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    def test_score_interface(self):
        pdb_file = "../test/data/9had_AB.pdb" # Needs to be A_B
        interface_scores = score_interface(pdb_file, binder_chain="B")

        print("TEST_SCORE_INTERFACE:", interface_scores)

        # These aren't 'ideal' values, just what was observed when the test was committed.
        self.assertAlmostEqual(interface_scores['binder_score'], 107.08)
        self.assertAlmostEqual(interface_scores['surface_hydrophobicity'], 0.16)
        self.assertAlmostEqual(interface_scores['interface_sc'], 0.69)

        # This is reporting 0.48 with PyRosetta4.conda.ubuntu.cxx11thread.serialization.Ubuntu.python312.Release 2025.17+release.356248d2035a0749e09a4a79479678a8f54e7220
        # 0.52 with PyRosetta4.conda.ubuntu-20.04.cxx11thread.serialization.Ubuntu.python314.Release 2026.33+release.97cfde1af8f45eee9625c0fd72c6b7194bf31d24
        # but 0.55 on the test server (Rosetta devel 2026.33.post.dev+2.HEAD.5c10c92fe4852989385dd261fa3e516a762bc395)
        self.assertGreater(interface_scores['interface_packstat'], 0.4)
        self.assertLess(interface_scores['interface_packstat'], 0.6)

        self.assertAlmostEqual(interface_scores['interface_dG'], 113.75)
        self.assertAlmostEqual(interface_scores['interface_dSASA'], 1853.22)
        self.assertAlmostEqual(interface_scores['interface_dG_SASA_ratio'], 6.14)
        self.assertAlmostEqual(interface_scores['interface_fraction'], 31.9)
        self.assertAlmostEqual(interface_scores['interface_nres'], 50)
        self.assertAlmostEqual(interface_scores['interface_interface_hbonds'], 5)
        self.assertAlmostEqual(interface_scores['interface_hbond_percentage'], 10.0)
        self.assertAlmostEqual(interface_scores['interface_delta_unsat_hbonds'], 4.0)
        self.assertAlmostEqual(interface_scores['interface_delta_unsat_hbonds_percentage'], 8.0)


    def test_align(self):
        reference_pdb = "../test/data/9had_AB.pdb"
        align_pdb = os.path.join(self.workdir.name, "align_test.pdb")
        shutil.copyfile("../test/data/9had_EF.pdb", align_pdb )
        reference_chain_id = "B"
        align_chain_id = "F"
        reference_chain_id_to = "A"
        align_chain_id_to = "E"

        pre_rmsd = unaligned_rmsd(reference_pdb, align_pdb, reference_chain_id, align_chain_id)

        align_pdbs(reference_pdb, align_pdb, reference_chain_id_to, align_chain_id_to)
        # This overwrites the align_pdb file contents.

        post_rmsd = unaligned_rmsd(reference_pdb, align_pdb, reference_chain_id, align_chain_id)

        print("ALIGN PRE:", pre_rmsd, "POST:", post_rmsd)

        self.assertNotEqual(pre_rmsd, post_rmsd)

        self.assertNotEqual(pre_rmsd, 0)
        self.assertNotEqual(post_rmsd, 0)
        self.assertGreater( pre_rmsd, post_rmsd )

        self.assertGreater( pre_rmsd, 40 ) # 47.96 when checked in
        self.assertLess( post_rmsd, 5 ) # 2.78 when checked in


    def test_pr_relax(self):
        pdb_file = "../test/data/9had_AB.pdb"

        scorefxn = pyrosetta.get_fa_scorefxn()

        relaxed_pdb_path = os.path.join(self.workdir.name, "relaxed.pdb")

        pose = rosetta.core.import_pose.pose_from_file(pdb_file)
        pre_score = scorefxn(pose)

        pr_relax(pdb_file, relaxed_pdb_path)

        pose = rosetta.core.import_pose.pose_from_file(relaxed_pdb_path)
        post_score = scorefxn(pose)

        print("RELAX PRE:", pre_score, "POST:", post_score)

        self.assertLess( post_score,  pre_score )

        self.assertGreater( pre_score , 0 )
        self.assertLess( post_score, 0 )

        self.assertGreater( pre_score , 200 ) # 234.70378992853523 when committed
        self.assertLess( post_score, -550 ) # -580.1071282502197 when committed


if __name__ == "__main__":
    unittest.main()
