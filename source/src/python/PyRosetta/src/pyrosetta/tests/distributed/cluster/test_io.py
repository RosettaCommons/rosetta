"""
PyRosettaCluster smoke tests using the `unittest` framework.
"""
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import glob
import json
import os
import pyrosetta.distributed
import tempfile
import unittest

from pyrosetta.distributed.cluster import run
from pyrosetta.distributed.cluster.io import secure_read_pickle


class IOTest(unittest.TestCase):
    _my_string_value = "foo"
    _my_real_value = 12.34567890123456789
    _my_pose_value = "DATA"
    _my_complex_value = 4j

    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        cls.workdir = tempfile.TemporaryDirectory()
        cls.decoy_dir_name = "decoys"
        cls.simulation_name = "TestIO"
        cls.instance_kwargs = dict(
            tasks=IOTest.create_tasks,
            protocols=IOTest.my_pyrosetta_protocol,
            input_packed_pose=None,
            seeds=None,
            decoy_ids=None,
            client=None,
            scheduler=None,
            scratch_dir=cls.workdir.name,
            cores=None,
            processes=None,
            memory=None,
            min_workers=1,
            max_workers=1,
            nstruct=1,
            dashboard_address=None,
            compression=True,
            logging_level="INFO",
            project_name="PyRosettaCluster_Tests",
            simulation_name=cls.simulation_name,
            environment=None,
            decoy_dir_name=cls.decoy_dir_name,
            logs_dir_name="logs",
            ignore_errors=False,
            timeout=0.1,
            max_delay_time=0.5,
            sha1=None,
            dry_run=False,
            save_all=False,
            system_info=None,
            pyrosetta_build=None,
        )

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    @staticmethod
    def my_pyrosetta_protocol(packed_pose, **kwargs):
        import pyrosetta
        import pyrosetta.distributed.io as io

        packed_pose = io.pose_from_sequence(kwargs["seq"])
        packed_pose = packed_pose.update_scores(
            my_string_score=IOTest._my_string_value,
            my_real_score=IOTest._my_real_value,
            my_pose_score=pyrosetta.pose_from_sequence(IOTest._my_pose_value),
            my_complex_score=IOTest._my_complex_value,
        )

        return packed_pose

    @staticmethod
    def create_tasks():
        tasks = []
        for i in range(1, 5):
            task = {
                "extra_options": "-ex1 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
                "seq": "TEST" * i,
            }
            tasks.append(task)
        return tasks

    def test_io(self):
        """Smoke test for basic PyRosettaCluster I/O."""
        output_decoy_types = [".pdb", ".b64_pose"]
        output_scorefile_types = [".json", ".gz", ".bz2", ".xz", ".tar.gz", ".tar.bz2", ".tar.xz", ".zip"]
        for compressed in (True, False):
            for simulation_records_in_scorefile in (True, False):
                for scorefile_name in (None, "test_my_scorefile.json"):
                    output_path = os.path.join(
                        self.workdir.name,
                        "outputs_{0}_{1}_{2}".format(compressed, simulation_records_in_scorefile, scorefile_name)
                    )
                    instance_kwargs = {
                        **self.instance_kwargs,
                        "compressed": compressed,
                        "simulation_records_in_scorefile": simulation_records_in_scorefile,
                        "scorefile_name": scorefile_name,
                        "output_decoy_types": output_decoy_types,
                        "output_scorefile_types": output_scorefile_types,
                        "output_path": output_path,
                    }
                    if "pandas" not in pyrosetta.secure_unpickle.get_secure_packages():
                        with self.assertRaises(AssertionError):  # output_scorefile_types=[".gz", ...] requires 'pandas' as a secure package
                            run(**instance_kwargs)
                    pyrosetta.secure_unpickle.add_secure_package("pandas")
                    run(**instance_kwargs)
                    # Test decoy outputs
                    _n_tasks = len(IOTest.create_tasks())
                    for output_decoy_type in output_decoy_types:
                        output_files = glob.glob(
                            os.path.join(
                                output_path,
                                self.decoy_dir_name,
                                "*",
                                f"*{output_decoy_type}.bz2" if compressed else f"*{output_decoy_type}",
                            )
                        )
                        self.assertEqual(_n_tasks, len(output_files), msg=f"Incorrect number of output files for type: {output_decoy_type}")
                    # Test scorefile outputs
                    for output_scorefile_type in output_scorefile_types:
                        _scorefile_name = "scores.json" if scorefile_name is None else scorefile_name
                        scorefile = os.path.join(output_path, _scorefile_name.replace(".json", output_scorefile_type))
                        self.assertTrue(os.path.isfile(scorefile), msg=f"Scorefile was not saved: {scorefile}")
                        if output_scorefile_type == ".json":
                            with open(scorefile, "r") as f:
                                entries = list(map(json.loads, f))
                            for entry in entries:
                                if simulation_records_in_scorefile:
                                    self.assertIn("instance", entry.keys())
                                    self.assertIn("metadata", entry.keys())
                                    self.assertIn("scores", entry.keys())
                                    self.assertIn("my_string_score", entry["scores"])
                                    self.assertEqual(entry["scores"]["my_string_score"], IOTest._my_string_value)
                                    self.assertIn("my_real_score", entry["scores"])
                                    self.assertEqual(entry["scores"]["my_real_score"], IOTest._my_real_value)
                                    self.assertNotIn("my_pose_score", entry["scores"])
                                    self.assertNotIn("my_complex_score", entry["scores"])
                                else:
                                    self.assertEqual(len(entry), 1)
                                    output_file, scores = next(iter(entry.items()))
                                    self.assertTrue(os.path.basename(output_file).startswith(self.simulation_name))
                                    self.assertIn("my_string_score", scores)
                                    self.assertEqual(scores["my_string_score"], IOTest._my_string_value)
                                    self.assertIn("my_real_score", scores)
                                    self.assertEqual(scores["my_real_score"], IOTest._my_real_value)
                                    self.assertNotIn("my_pose_score", scores)
                                    self.assertNotIn("my_complex_score", scores)
                        else:
                            df = secure_read_pickle(scorefile, compression="infer")
                            if simulation_records_in_scorefile:
                                self.assertIn("instance", df.columns)
                                self.assertIn("metadata", df.columns)
                                self.assertIn("scores", df.columns)
                                scores = df["scores"]
                                for index in scores.index:
                                    self.assertIn("my_string_score", scores.loc[index].keys())
                                    self.assertEqual(scores.loc[index]["my_string_score"], IOTest._my_string_value)
                                    self.assertIn("my_real_score", scores.loc[index].keys())
                                    self.assertEqual(scores.loc[index]["my_real_score"], IOTest._my_real_value)
                                    self.assertIn("my_pose_score", scores.loc[index].keys())
                                    self.assertEqual(scores.loc[index]["my_pose_score"].sequence(), IOTest._my_pose_value)
                                    self.assertIn("my_complex_score", scores.loc[index].keys())
                                    self.assertEqual(scores.loc[index]["my_complex_score"], IOTest._my_complex_value)
                            else:
                                for index in df.index:
                                    self.assertIn("my_string_score", df.columns)
                                    self.assertEqual(df.at[index, "my_string_score"], IOTest._my_string_value)
                                    self.assertIn("my_real_score", df.columns)
                                    self.assertEqual(df.at[index, "my_real_score"], IOTest._my_real_value)
                                    self.assertIn("my_pose_score", df.columns)
                                    self.assertEqual(df.at[index, "my_pose_score"].sequence(), IOTest._my_pose_value)
                                    self.assertIn("my_complex_score", df.columns)
                                    self.assertEqual(df.at[index, "my_complex_score"], IOTest._my_complex_value)
