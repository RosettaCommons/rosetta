# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import glob
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose
import pyrosetta.distributed.tasks.score as score
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
import tempfile
import unittest 


class SmokeTestDistributed(unittest.TestCase):

    min_rs = """
        <ROSETTASCRIPTS>
            <SCOREFXNS>
            </SCOREFXNS>

            <RESIDUE_SELECTORS>
            </RESIDUE_SELECTORS>

            <TASKOPERATIONS>
            </TASKOPERATIONS>

            <FILTERS>
            </FILTERS>

            <MOVERS>
            </MOVERS>

            <PROTOCOLS>
            </PROTOCOLS>

            <OUTPUT />
        </ROSETTASCRIPTS>
    """

    def test_score_smoke_test(self):
        """RosettaScripts 'null' tasks just score with default score function.

        A high-level smoke test of the distributed namespace. Inits a
        packed-pose object via a call through the io layer, passes this through
        the RosettaScripts layer, then through the scoring layers. Covers pose
        serializiation/deserialization, score value extraction and rosetta
        scripts parser access.

        Which is to say, turn on the power and look for magic smoke.
        """
        flags = """
            -ignore_unrecognized_res  1
               -ex4 # Test comment 1
             -out:level     300   ### Test comment 2
            """
        pyrosetta.distributed.init(flags)

        flags_dict = {
            "-ignore_unrecognized_res": "1",  # Test comment 1
            "-out:level": "500",  # Test comment 2
        }
        self.assertEqual(
            pyrosetta.distributed._normflags(flags_dict),
            "-ignore_unrecognized_res 1 -out:level 500",
        )

        score_task = score.ScorePoseTask()
        rs_task = rosetta_scripts.SingleoutputRosettaScriptsTask(self.min_rs)

        score_result = score_task(io.pose_from_sequence("TEST")).scores
        rs_result = rs_task(io.pose_from_sequence("TEST")).scores

        self.assertAlmostEqual(
            score_result["total_score"], rs_result["total_score"]
        )

        score_task = io.create_score_function("score12")
        self.assertIsInstance(score_task, score.ScorePoseTask)

        score_task = io.get_fa_scorefxn()
        self.assertIsInstance(score_task, score.ScorePoseTask)

        score_task = io.get_score_function()
        self.assertIsInstance(score_task, score.ScorePoseTask)

    def test_update_score(self):
        """PackedPose.update_score returns an updated copy.

        PackedPose.update_score performs an copy-update of the pack, not an
        inplace modification. New score values are applied from kwargs and args,
        with kwargs and later-arg masking duplicate values.
        """

        work_pose = io.pose_from_sequence("TEST")

        self.assertDictEqual(work_pose.scores, dict())

        # Test merge and masking, just args
        work_updated = work_pose.update_scores(
            {"arg": 1.0, "dupe": "foo"},
            dict(arg2=2.0, dupe="bar")
        )

        self.assertDictEqual(work_pose.scores, dict())
        self.assertDictEqual(work_updated.scores, dict(arg=1.0, arg2=2.0, dupe="bar"))

        # Test merge and masking, args and kwargs
        work_updated = work_pose.update_scores(
            {"arg": 1.0, "dupe": "foo"},
            dict(arg2=2.0, dupe="bar"),
            kwarg="yes",
            dupe="bat"
        )

        self.assertDictEqual(work_pose.scores, dict())
        self.assertDictEqual(work_updated.scores, dict(arg=1.0, arg2=2.0, kwarg="yes", dupe="bat"))

        # Test just kwargs
        work_updated = work_pose.update_scores(
            kwarg="yes",
            dupe="bat"
        )

        self.assertDictEqual(work_pose.scores, dict())
        self.assertDictEqual(work_updated.scores, dict(kwarg="yes", dupe="bat"))

        # Test no args
        work_updated = work_pose.update_scores()

        self.assertDictEqual(work_pose.scores, dict())
        self.assertDictEqual(work_updated.scores, dict())

    def test_silent_io(self):

        with tempfile.TemporaryDirectory() as workdir:

            # Tests if single pose is written to silent file and returned unchanged.
            test_packed_pose = io.pose_from_sequence("TESTTESTTEST")
            test_pose = packed_pose.to_pose(test_packed_pose)
            tmp_file = os.path.join(workdir, "temp.silent")
            io.to_silent(test_packed_pose, tmp_file)
            returned_poses = [packed_pose.to_pose(p) for p in io.poses_from_silent(tmp_file)]
            
            # Tests that the amino acid sequences does not change after 
            # being written to silent file and read back in.
            # Cannot just assert that the two poses are equal, 
            # because writing them to the silent file adds score lines. 
            self.assertEqual(
                test_pose.sequence(), 
                returned_poses[0].sequence(), 
                msg="Single sequence recovery failed.")

            # Tests that the positions of atoms are almost identical after
            # being written to silent file and read back in.
            # Rmsd will not quite be equal because the output silent file 
            # truncates xyz coordinates to the third decimal place.
            self.assertAlmostEqual(
                0., 
                pyrosetta.rosetta.core.scoring.all_atom_rmsd(test_pose, returned_poses[0]), 
                places=3,
                msg="Single position recovery failed.")

            # Test if a list of poses can be written to a silent file 
            # and returned unchanged.
            test_packed_poses = [io.pose_from_sequence("TEST" * i) for i in range(1, 5)]
            test_poses = [packed_pose.to_pose(p) for p in test_packed_poses]
            tmp_file = os.path.join(workdir, "temp_list.silent")
            io.to_silent(test_packed_poses, tmp_file)
            returned_poses = [packed_pose.to_pose(p) for p in io.poses_from_silent(tmp_file)]
            for i in range(len(test_poses)):
                self.assertEqual(
                    test_poses[i].sequence(), 
                    returned_poses[i].sequence(), 
                    msg="List sequence recovery failed.")
                self.assertAlmostEqual(
                    0.,
                    pyrosetta.rosetta.core.scoring.all_atom_rmsd(test_poses[i],returned_poses[i]), 
                    places=3,
                    msg="List position recovery failed.")

    def test_clone(self):
        # Test PackedPose cloning
        packed_pose_1 = io.pose_from_sequence("DNA")
        packed_pose_1.scores["key"] = 37.0
        packed_pose_2 = packed_pose_1.clone()
        self.assertNotEqual(
            id(packed_pose_1),
            id(packed_pose_2),
            msg="Clone failed.",
        )
        self.assertEqual(
            packed_pose_1.scores["key"],
            packed_pose_2.scores["key"],
            msg="Scores differ after clone.",
        )
        self.assertNotEqual(
            id(packed_pose_1.scores["key"]),
            id(packed_pose_2.scores["key"]),
            msg="Score memory address is identical after clone.",
        )

    def roundtrip(self, func, ext, input_packed_pose, workdir, scorefxn):
        """Used in `test_packed_pose_io` testing framework."""
        # Dump PackedPose to disk
        out_file = os.path.join(workdir, f"tmp.{ext}")
        if ext == "scored.pdb":
            func(input_packed_pose.pose, out_file, "score12")
            func(input_packed_pose, out_file, score.ScorePoseTask(weights="ref2015_cart"))
            func(input_packed_pose, out_file, scorefxn)
        else:
            func(input_packed_pose, out_file)
        # Load PackedPose from disk
        if ext in ("pdb.bz2", "bz2", "pdb.gz", "gz", "pdb.xz", "xz"):
            with self.assertRaises(FileNotFoundError):
                io.pose_from_pdb(123)
            with self.assertRaises(FileNotFoundError):
                io.pose_from_pdb(os.path.join(workdir, f"nonexistent_file.{ext}"))
            self.assertIsNone(io.pose_from_pdb(None))
            output_packed_pose = io.pose_from_pdb(out_file)
        else:
            output_packed_pose = io.pose_from_file(out_file)
        # Test annotated sequence recovery
        self.assertEqual(
            input_packed_pose.pose.annotated_sequence(),
            output_packed_pose.pose.annotated_sequence(),
            msg=f"Sequence recovery failed for extension '{ext}'."
        )
        # Test coordiante recovery
        places = 32 if ext in ("b64", "base64", "B64", "pose", "pickle", "pickled_pose") else 2
        for res in range(1, input_packed_pose.pose.size() + 1):
            for atom in range(1, input_packed_pose.pose.residue(res).natoms() + 1):
                atom_input = input_packed_pose.pose.residue(res).xyz(atom)
                atom_output = output_packed_pose.pose.residue(res).xyz(atom)
                for axis in "xyz":
                    self.assertAlmostEqual(
                        getattr(atom_input, axis),
                        getattr(atom_output, axis),
                        places=places,
                        msg=f"Coordinate recovery failed for extension '{ext}'."
                    )
        # Test score recovery
        self.assertTrue(dict(input_packed_pose.pose.scores), msg=f"Pose scores dictionary is empty for extension '{ext}'.")
        if ext in ("pdb", "scored.pdb", "cif", "mmcif", "mmtf", "pdb.bz2", "bz2", "pdb.gz", "gz", "pdb.xz", "xz"):
            self.assertFalse(dict(output_packed_pose.pose.scores), msg=f"Pose scores dictionary has items for extension '{ext}'.")
            self.assertNotEqual(
                input_packed_pose.scores,
                output_packed_pose.scores,
                msg=f"PackedPose scores dictionaries differ for extension '{ext}'."
            )
        else: # base64-encoded and pickle-encoded files save the `pose.scores` dictionary
            self.assertTrue(output_packed_pose.scores, msg=f"PackedPose scores dictionary is empty for extension '{ext}'.")
            self.assertTrue(dict(output_packed_pose.pose.scores), msg=f"Pose scores dictionary is empty for extension '{ext}'.")
            self.assertDictEqual(
                input_packed_pose.scores,
                output_packed_pose.scores,
                msg=f"Pose scores dictionaries differ for extension '{ext}'."
            )
            self.assertDictEqual(
                dict(input_packed_pose.pose.scores),
                dict(output_packed_pose.pose.scores),
                msg=f"Pose scores dictionaries differ for extension '{ext}'."
            )

    def test_packed_pose_io(self):
        # Test PackedPose I/O
        ext_func_dict = {
            "pdb": io.dump_pdb,
            "scored.pdb": io.dump_scored_pdb,
            "cif": io.dump_cif,
            "mmcif": io.dump_cif,
            "mmtf": io.dump_mmtf,
            "pdb.bz2": io.dump_pdb_bz2,
            "bz2": io.dump_pdb_bz2,
            "pdb.gz": io.dump_pdb_gz,
            "gz": io.dump_pdb_gz,
            "pdb.xz": io.dump_pdb_xz,
            "xz": io.dump_pdb_xz,
            "base64": io.dump_base64,
            "b64": io.dump_base64,
            "B64": io.dump_base64,
            "pose": io.dump_base64,
            "pickle": io.dump_pickle,
            "pickled_pose": io.dump_pickle,
        }
        input_packed_pose = io.pose_from_sequence("CATALYST/X[ATP]")
        scorefxn = pyrosetta.create_score_function("ref2015")
        input_pose = input_packed_pose.pose
        scorefxn(input_pose)
        input_packed_pose = io.to_packed(input_pose)
        with tempfile.TemporaryDirectory() as workdir:
            for ext, func in ext_func_dict.items():
                self.roundtrip(func, ext, input_packed_pose.clone(), workdir, scorefxn)

    def test_poses_from_sequences_io(self):
        # Test `io.poses_from_sequences`
        seqs = ["DESIGN" * i for i in range(3, 8)]
        in_packed_poses = list(io.poses_from_sequences(seqs))
        for packed_pose, seq in zip(in_packed_poses, seqs):
            self.assertIsInstance(
                packed_pose,
                pyrosetta.distributed.packed_pose.core.PackedPose,
                msg="Instance is not a `PackedPose` object."
            )
            self.assertEqual(packed_pose.pose.sequence(), seq, msg="Sequences diverge.")

    def test_multimodel_pdb_io(self):
         # Test `poses_from_multimodel_pdb`
        seqs = ["PHENYLALANINE" * i for i in range(3, 8)]
        in_packed_poses = list(io.poses_from_sequences(seqs))
        with tempfile.TemporaryDirectory() as workdir:
            pdb_file = os.path.join(workdir, "tmp.pdb")
            io.dump_multimodel_pdb(in_packed_poses, pdb_file)
            out_packed_poses = list(io.poses_from_multimodel_pdb(pdb_file))
        self.assertEqual(
            len(out_packed_poses),
            len(in_packed_poses),
            msg="Number of `PackedPose` objects recovery failed."
        )
        for in_packed_pose, out_packed_pose in zip(in_packed_poses, out_packed_poses):
            self.assertIsInstance(
                in_packed_pose,
                pyrosetta.distributed.packed_pose.core.PackedPose,
                msg="Instance is not a `PackedPose` object."
            )
            self.assertIsInstance(
                out_packed_pose,
                pyrosetta.distributed.packed_pose.core.PackedPose,
                msg="Instance is not a `PackedPose` object."
            )
            self.assertEqual(
                out_packed_pose.pose.annotated_sequence(),
                in_packed_pose.pose.annotated_sequence(),
                msg="Annotated sequence recovery failed."
            )

    def test_poses_from_files_io(self):
        seqs = ["/".join(["PHI/PSI"] * i) for i in range(1, 5)]
        in_packed_poses = io.poses_from_sequences(seqs)
        # Test `io.poses_from_files`
        with tempfile.TemporaryDirectory() as workdir:
            for i, test_packed_pose in enumerate(in_packed_poses, start=1):
                io.dump_pdb(test_packed_pose, os.path.join(workdir, "%s.pdb" % i))
            out_packed_poses = io.poses_from_files(glob.glob(os.path.join(workdir, "*.pdb")))
        for in_packed_pose, out_packed_pose in zip(in_packed_poses, out_packed_poses):
            self.assertIsInstance(
                in_packed_pose,
                pyrosetta.distributed.packed_pose.core.PackedPose,
                msg="Instance is not a `PackedPose` object."
            )
            self.assertIsInstance(
                out_packed_pose,
                pyrosetta.distributed.packed_pose.core.PackedPose,
                msg="Instance is not a `PackedPose` object."
            )
            self.assertEqual(
                out_packed_pose.pose.annotated_sequence(),
                in_packed_pose.pose.annotated_sequence(),
                msg="Annotated sequence recovery failed."
            )
