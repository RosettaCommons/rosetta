# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

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

        score_task = score.ScorePoseTask()
        rs_task = rosetta_scripts.SingleoutputRosettaScriptsTask(self.min_rs)

        score_result = score_task(io.pose_from_sequence("TEST")).scores
        rs_result = rs_task(io.pose_from_sequence("TEST")).scores

        self.assertAlmostEqual(
            score_result["total_score"], rs_result["total_score"]
        )

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
