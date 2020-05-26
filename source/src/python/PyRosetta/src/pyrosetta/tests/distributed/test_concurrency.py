# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import numpy
import time
import threading
import unittest

import pyrosetta.distributed
import pyrosetta.distributed.io as io
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
import pyrosetta.distributed.tasks.score as score

pyrosetta.distributed.init(options="-run:constant_seed 1", set_logging_handler=None)
# TODO: In multithreaded builds the python-based logging handler appears to fail using concurrent.futures...

class TestConcurrentScripts(unittest.TestCase):

    def test_concurrent_on_task(self):

        protocol = rosetta_scripts.SingleoutputRosettaScriptsTask("""
        <ROSETTASCRIPTS>
            <TASKOPERATIONS>
                #We'll only allow two packing threads within the three launched by this test (for a total of 6), to
                #avoid a thread explosion:
                <RestrictInteractionGraphThreadsOperation name="restrict_threads" thread_limit="4" />
            </TASKOPERATIONS>
            <MOVERS>
                <FastRelax name="frlx" repeats="1" task_operations="restrict_threads" />
            </MOVERS>
            <PROTOCOLS>
                <Add mover_name="frlx"/>
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
                <TASKOPERATIONS>
                    #We'll only allow two packing threads within the three launched by this test (for a total of 6), to
                    #avoid a thread explosion:
                    <RestrictInteractionGraphThreadsOperation name="restrict_threads" thread_limit="4" />
                </TASKOPERATIONS>
                <MOVERS>
                    <FastRelax name="frlx" repeats="1" task_operations="restrict_threads" />
                </MOVERS>
                <PROTOCOLS>
                    <Add mover_name="frlx"/>
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
