# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import argparse
import dask.distributed
import logging
import tempfile
import unittest

import pyrosetta.distributed.io as io
import pyrosetta.distributed.packed_pose as packed_pose
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts


class TestDaskDistribution(unittest.TestCase):
    
    _dask_scheduler = None
    
    with tempfile.TemporaryDirectory() as workdir:
        
        def setUp(self, local_dir=workdir):
            if not self._dask_scheduler:
                self.local_cluster = dask.distributed.LocalCluster(
                    n_workers=2, threads_per_worker=2, diagnostics_port=None, local_dir=local_dir
                )
                cluster = self.local_cluster
            else:
                self.local_cluster = None
                cluster = self._dask_scheduler

            self.client = dask.distributed.Client(cluster)

        def tearDown(self):
            self.client.close()
            
            if self.local_cluster:
                self.local_cluster.close()

        def test_rosetta_scripts(self):
            test_protocol = """
            <ROSETTASCRIPTS>
                <TASKOPERATIONS>
                    <RestrictToRepacking name="repack"/>
                </TASKOPERATIONS>
                <MOVERS>
                    <PackRotamersMover name="pack" task_operations="repack"/>
                </MOVERS>
                <PROTOCOLS>
                    <Add mover="pack"/>
                </PROTOCOLS>
            </ROSETTASCRIPTS>
            """

            test_pose = io.pose_from_sequence("TEST")
            test_task = rosetta_scripts.SingleoutputRosettaScriptsTask(test_protocol)

            logging.info("dask client: %s", self.client)
            task = self.client.submit(test_task, test_pose)
            result = task.result()
            self.assertEqual(
                packed_pose.to_pose(result).sequence(),
                packed_pose.to_pose(test_pose).sequence()
            )

        def test_basic(self):
            logging.info("dask client: %s", self.client)
            task = self.client.submit(lambda a: a + 1, 1)
            self.assertEqual(task.result(), 2)


if __name__ == "__main__":
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s.%(msecs).03d %(name)s %(message)s",
        datefmt='%Y-%m-%dT%H:%M:%S'
    )

    parser = argparse.ArgumentParser(
        description="Run initial pyrosetta.distributed smoke test over given scheduler."
    )
    parser.add_argument("scheduler", type=str, nargs="?", help="Target scheduler endpoint for test.")
    args = parser.parse_args()

    TestDaskDistribution._dask_scheduler = args.scheduler

    test_suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestDaskDistribution)
    runner = unittest.TextTestRunner()
    runner.run(test_suite)
