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
