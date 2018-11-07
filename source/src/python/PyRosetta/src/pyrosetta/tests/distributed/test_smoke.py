import unittest 

import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import pyrosetta.distributed.tasks.score as score
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts


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
