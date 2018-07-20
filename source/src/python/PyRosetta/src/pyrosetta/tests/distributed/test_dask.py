import logging

import argparse
import unittest

import dask.distributed

import pyrosetta.distributed.io as io
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts

class TestDaskDistribution(unittest.TestCase):
    _dask_scheduler = None

    def setUp(self):
        if not self._dask_scheduler:
            self.local_cluster = dask.distributed.LocalCluster(n_workers=2, threads_per_worker=2)
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
            <MOVERS>
                <PackRotamersMover name="pack"/>
            </MOVERS>
            <PROTOCOLS>
                <Add mover="pack" />
            </PROTOCOLS>
        </ROSETTASCRIPTS>
        """

        test_pose = io.pose_from_sequence("TESTTESTTEST")
        test_task = rosetta_scripts.SingleoutputRosettaScriptsTask(test_protocol)

        logging.info("dask client: %s", self.client)
        task = self.client.submit(test_task, test_pose)
        result = task.result()

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

    parser = argparse.ArgumentParser(description='Run initial pyrosetta.distributed smoke test over given scheduler.')
    parser.add_argument('scheduler', type=str, nargs='?', help='Target scheduler endpoint for test.')
    args = parser.parse_args()

    TestDaskDistribution._dask_scheduler = args.scheduler

    test_suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestDaskDistribution)
    runner = unittest.TextTestRunner()
    runner.run(test_suite)
