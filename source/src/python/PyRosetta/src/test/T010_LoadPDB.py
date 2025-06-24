#!/usr/bin/env python
# :noTabs=true:
# -*- coding: utf-8 -*-

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   LoadPDB.py
## @brief  Serve as example on how to import rosetta, load pdb and as script that create database
## @brief        binaries on windows.
## @author Sergey Lyskov
## @author Jason C. Klima

from __future__ import print_function

import os
import pyrosetta
import pyrosetta.rosetta as rosetta
import tempfile
import unittest

class LoadPDBTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
        print( pyrosetta.version() )
        cls.pose = pyrosetta.pose_from_sequence("DSEEKFLRRIGRFGYGYGPYE")
        cls.scorefxn = rosetta.core.scoring.get_score_function()
        cls.workdir = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    def test_load_pdb(self):
        os.chdir('.test.output')
        pose = rosetta.core.import_pose.pose_from_file("../test/data/test_in.pdb")

        scorefxn = rosetta.core.scoring.get_score_function()
        scorefxn(pose)


        pose2 = pyrosetta.pose_from_sequence("ARNDCEQGHILKMFPSTWYV", 'fa_standard')

        scorefxn = rosetta.core.scoring.get_score_function()
        scorefxn(pose2)


        pose3 = pyrosetta.pose_from_sequence("DSEEKFLRRIGRFGYGYGPYE",'centroid')

        # Creating standard centroid score function and scoring
        scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function('score3')
        scorefxn(pose3)

        pose_fs = pyrosetta.pose_from_sequence("DSEEKFLRRIGRFGYGYGPYE")
        pose_fs.delete_polymer_residue(2)  # Testing that attached PDB info have right size...


        poses = pyrosetta.poses_from_silent('../test/data/test_in.silent')
        for pose in poses:
            scorefxn = rosetta.core.scoring.get_score_function()
            scorefxn(pose)

        import pyrosetta.toolbox

        # commenting this out for now to avoid release failures during debug builds when network is out
        # pyrosetta.toolbox.pose_from_rcsb('1brs')

    def test_roundtrip_pdb(self):
        for ext in (".pdb", ".pdb.bz2", ".bz2", ".pdb.gz", ".gz", ".pdb.xz", ".xz"):
            pdb_file = os.path.join(self.workdir.name, "tmp_pdb{0}".format(ext))
            pyrosetta.dump_pdb(self.pose, pdb_file)
            pose_out = pyrosetta.pose_from_file(pdb_file)
            self.assertEqual(pose_out.size(), self.pose.size())

    def test_roundtrip_scored_pdb(self):
        for ext in (".pdb", ".pdb.bz2", ".bz2", ".pdb.gz", ".gz", ".pdb.xz", ".xz"):
            pdb_file = os.path.join(self.workdir.name, "tmp_scored_pdb{0}".format(ext))
            pyrosetta.dump_scored_pdb(self.pose, pdb_file, self.scorefxn)
            pose_out = pyrosetta.pose_from_file(pdb_file)
            self.assertEqual(pose_out.size(), self.pose.size())

    def test_roundtrip_file(self):
        for ext in (".cif", ".mmcif", ".mmtf"):
            pdb_file = os.path.join(self.workdir.name, "tmp_file{0}".format(ext))
            pyrosetta.dump_file(self.pose, pdb_file)
            pose_out = pyrosetta.pose_from_file(pdb_file)
            self.assertEqual(pose_out.size(), self.pose.size())

    def test_pdbstring_io(self):
        self.assertIn("ATOM", pyrosetta.io.to_pdbstring(self.pose))
        pose_out = pyrosetta.io.pose_from_pdbstring(pyrosetta.io.to_pdbstring(self.pose))
        self.assertEqual(pose_out.size(), self.pose.size())

if __name__ == "__main__":
    unittest.main()
