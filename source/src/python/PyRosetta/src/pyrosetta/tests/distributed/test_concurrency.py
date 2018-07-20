import unittest
import threading
import time
import numpy

import pyrosetta.distributed.io as io
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
import pyrosetta.distributed.tasks.score as score


class TestConcurrentScripts(unittest.TestCase):
    def test_concurrent_on_task(self):

        protocol = rosetta_scripts.SingleoutputRosettaScriptsTask("""
        <ROSETTASCRIPTS>

        <MOVERS>
            <FastRelax name="score" repeats="1" />
        </MOVERS>

        <PROTOCOLS>
            <Add mover_name="score"/>
        </PROTOCOLS>

        </ROSETTASCRIPTS>
        """)

        test_pose = io.pose_from_sequence("TEST")

        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as p:
            result = list(p.map(protocol, [test_pose] * 3))

    def test_concurrent_multi_task(self):

        def run_task(seq):
            test_pose = io.pose_from_sequence(seq)
            protocol = rosetta_scripts.SingleoutputRosettaScriptsTask("""
            <ROSETTASCRIPTS>

            <MOVERS>
                <FastRelax name="score" repeats="1" />
            </MOVERS>

            <PROTOCOLS>
                <Add mover_name="score"/>
            </PROTOCOLS>

            </ROSETTASCRIPTS>
            """)

            return protocol(test_pose)


        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as p:
            result = list(p.map(run_task, ["TEST"] * 3))

    def test_concurrent_score(self):

        test_pose = io.pose_from_sequence("TEST" * 10)

        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as p:
            result = p.map(score.ScorePoseTask(), [test_pose] * 3)
