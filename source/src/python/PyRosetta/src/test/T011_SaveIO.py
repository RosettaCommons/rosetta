#!/usr/bin/env python
# :noTabs=true:
# -*- coding: utf-8 -*-

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   T011_SavePDB.py
## @brief  Tests for the saving functionality of PyRosetta
## @author Rocco Moretti

from __future__ import print_function

import os
import pyrosetta
import pyrosetta.rosetta as rosetta
import tempfile
import unittest
import json

try:
    import lzma as xz
    _skip_xz = False
except ImportError:
    _skip_xz = True


class SavePDBTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
        print( pyrosetta.version() )
        cls.pose1 = rosetta.core.import_pose.pose_from_file("../test/data/test_in.pdb")
        cls.pose2 = pyrosetta.pose_from_sequence("ARNDCEQGHILKMFPSTWYV", 'fa_standard')
        cls.scorefxn = rosetta.core.scoring.get_score_function()
        cls.scorefxn(cls.pose1)
        cls.scorefxn(cls.pose2)
        cls.workdir = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    def test_save_silent(self):

        silent_filename = os.path.join(self.workdir.name,"output.silent")

        pyrosetta.poses_to_silent( [self.pose1, self.pose2], silent_filename )

        self.assertTrue( os.path.exists(silent_filename) )

        num_score = 0
        num_annotated = 0
        with open(silent_filename) as f:
            for line in f:
                if line.startswith("SCORE:"):
                    num_score += 1
                elif line.startswith("ANNOTATED_SEQUENCE"):
                    num_annotated += 1
                # Other checks?

        self.assertEqual(num_score, 3) # plus the header
        self.assertEqual(num_annotated, 2)

    def test_rountrip_silent(self):
        silent_filename = os.path.join(self.workdir.name,"roundtrip.silent")

        pyrosetta.poses_to_silent( [self.pose1, self.pose2], silent_filename )

        self.assertTrue( os.path.exists(silent_filename) )

        poses = list( pyrosetta.poses_from_silent(silent_filename) )

        self.assertEquals(len(poses), 2 )
        self.assertEquals(poses[0].size(), self.pose1.size())
        self.assertEquals(poses[1].size(), self.pose2.size())
        self.assertEquals(poses[0].annotated_sequence(), self.pose1.annotated_sequence())
        self.assertEquals(poses[1].annotated_sequence(), self.pose2.annotated_sequence())

    def test_scorefile(self):
        scorefile_name = os.path.join(self.workdir.name,"scorefile.sc")

        pyrosetta.poses_to_scorefile( [self.pose2, self.pose1], scorefile_name ) )

        self.assertTrue( os.path.exists(scorefile_name) )

        lengths = []
        n_score = 0
        with open(scorefile_name) as f:
            for line in f:
                lengths.append( len(line.split()) )
                if line.startswith("SCORE:"):
                    n_score += 1

        self.assertEqual(n_score, 3) # plus the header
        self.assertTrue( len(lengths) <= n_score+1 ) # Allow for spurious "SEQUENCE" line
        lenghts = lengths[-n_score:]
        self.assertEqual( min(lenghts), max(lengths) )

    def test_json_scorefile(self):
        scorefile_name = os.path.join(self.workdir.name,"scorefile.sc")

        pose_clone = self.pose1.clone()
        pyrosetta.rosetta.core.pose.setPoseExraScore(pose_clone, "extra_real", 3.14159 )
        pyrosetta.rosetta.core.pose.setPoseExraScore(pose_clone, "extra_string", "TAG_VALUE" )

        pyrosetta.poses_to_scorefile( pose_clone, scorefile_name )
        pyrosetta.poses_to_scorefile( [self.pose2, self.pose1], scorefile_name ), use_json=True )

        self.assertTrue( os.path.exists(scorefile_name) )

        with open(scorefile_name) as f:
            lines = f.readlines()

        self.assertEqual(len(lines),3)

        has_extra_real = 0
        has_extra_string = 0
        for line in lines:
            try:
                scores = json.loads(line)
            except:
                self.assertTrue(False)
                continue

            self.assertTrue( scores.has("total_score") )
            if scores.has("extra_real"):
                has_extra_real += 1
            if scores.has("extra_string"):
                has_extra_string += 1

        self.assertEqual(has_extra_real,1)
        self.assertEqual(has_extra_string,1)

    def test_get_scorefile_info(self):
        pose_clone = self.pose1.clone()
        pyrosetta.rosetta.core.pose.setPoseExraScore(pose_clone, "extra_real", 3.14159 )
        pyrosetta.rosetta.core.pose.setPoseExraScore(pose_clone, "extra_string", "TAG_VALUE" )

        info = pyrosetta.io.get_scorefile_info(pose_clone)

        self.assertTrue( info.has("total_score") )
        self.assertAlmoseEqual( info["extra_real"], 3.14159 )
        self.assertEqual( info["extra_string"], "TAG_VALUE" )

if __name__ == "__main__":
    unittest.main()
