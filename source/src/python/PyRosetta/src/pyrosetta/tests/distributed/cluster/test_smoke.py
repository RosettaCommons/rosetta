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
__email__ = "klima.jason@gmail.com"

import json
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import random
import tempfile
import unittest

from pyrosetta.distributed.cluster import PyRosettaCluster, produce, reserve_scores, run


class SmokeTest(unittest.TestCase):
    def test_smoke(self):
        """Smoke test for basic PyRosettaCluster usage."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )

        def create_tasks():
            for i in range(1, 5):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "LYELL" * i,
                }

        def my_pyrosetta_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            pose_from_kwargs = pyrosetta.io.pose_from_sequence(kwargs["seq"])
            self.assertNotEqual(pose_from_kwargs.sequence(), "")
            self.assertNotIn("TESTING", pose_from_kwargs.sequence())
            pose_from_kwargs.clear()
            self.assertEqual(pose_from_kwargs.sequence(), "")
            pose_from_args = io.to_pose(packed_pose)
            self.assertEqual(pose_from_args.sequence(), "TESTING")
            self.assertNotIn("LYELL", pose_from_args.sequence())

            return packed_pose

        with tempfile.TemporaryDirectory() as workdir:

            instance_kwargs = dict(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence("TESTING"),
                seeds=None,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=1,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name=None,
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=os.path.join(workdir, "outputs"),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=0.1,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )
            cluster = PyRosettaCluster(**instance_kwargs)
            cluster.distribute(my_pyrosetta_protocol,)

            instance_kwargs.update({"protocols": my_pyrosetta_protocol})
            produce(**instance_kwargs)
            run(**instance_kwargs)


class SmokeTestMulti(unittest.TestCase):
    def test_smoke_multi(self):
        """Smoke test for PyRosettaCluster usage with multiple protocols."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )

        def create_tasks():
            for i in range(1, 5):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "TEST" * i,
                }

        def my_first_protocol(packed_pose, **kwargs):
            import pyrosetta

            pose = pyrosetta.io.pose_from_sequence(kwargs["seq"])
            pyrosetta.rosetta.core.pose.setPoseExtraScore(
                pose, "test_setPoseExtraScore", 123
            )

            return [pose.clone() for _ in range(3)]

        @reserve_scores
        def my_second_protocol(packed_pose, **kwargs):
            import pyrosetta  # noqa
            import pyrosetta.distributed.io as io  # noqa

            self.assertDictContainsSubset(
                {"test_setPoseExtraScore": 123}, packed_pose.scores
            )
            packed_pose.scores.clear()
            self.assertDictEqual({}, packed_pose.scores)
            pose = io.to_pose(packed_pose)
            for _ in range(3):
                yield pose.clone()

        def my_third_protocol(packed_pose, **kwargs):
            self.assertDictContainsSubset(
                {"test_setPoseExtraScore": 123}, packed_pose.scores
            )
            return my_second_protocol(packed_pose, **kwargs)

        with tempfile.TemporaryDirectory() as workdir:
            PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence("TESTING"),
                seeds=None,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=1,
                dashboard_address=None,
                compressed=False,
                logging_level="DEBUG",
                scorefile_name=None,
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=os.path.join(workdir, "outputs"),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=0.1,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            ).distribute(
                protocols=(my_first_protocol, my_second_protocol, my_third_protocol),
            )

    def test_smoke_multi_from_instance(self):
        """Smoke test for PyRosettaCluster usage with multiple protocols and instances."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )

        def my_first_protocol(packed_pose, **kwargs):
            import pyrosetta

            pose = pyrosetta.io.pose_from_sequence(kwargs["seq"])
            pyrosetta.rosetta.core.pose.setPoseExtraScore(
                pose, "test_setPoseExtraScore", 123
            )
            self.assertEqual(kwargs["PyRosettaCluster_protocol_number"], 0)
            return [pose.clone() for _ in range(3)]

        @reserve_scores
        def my_second_protocol(packed_pose, **kwargs):
            import pyrosetta  # noqa
            import pyrosetta.distributed.io as io

            self.assertDictContainsSubset(
                {"test_setPoseExtraScore": 123}, packed_pose.scores
            )
            packed_pose.scores.clear()
            self.assertDictEqual({}, packed_pose.scores)
            self.assertIn(kwargs["PyRosettaCluster_protocol_number"], [1, 2])
            pose = io.to_pose(packed_pose)
            for _ in range(3):
                yield pose.clone()

        def my_third_protocol(packed_pose, **kwargs):
            self.assertDictContainsSubset(
                {"test_setPoseExtraScore": 123}, packed_pose.scores
            )
            self.assertEqual(kwargs["PyRosettaCluster_protocol_number"], 2)
            return my_second_protocol(packed_pose, **kwargs)

        with tempfile.TemporaryDirectory() as workdir:
            PyRosettaCluster(
                tasks={
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "TEST" * 7,
                },
                input_packed_pose=io.pose_from_sequence("LYELL"),
                seeds=[
                    random.randint(-int((2 ** 32) / 2), int((2 ** 32) / 2) - 1)
                    for _ in range(3)
                ],
                decoy_ids=[random.randint(0, 2) for _ in range(3)],
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=12,
                nstruct=1,
                dashboard_address=None,
                compressed=True,
                logging_level="WARNING",
                scorefile_name=None,
                project_name="PyRosettaCluster_Testing",
                simulation_name=None,
                environment=None,
                output_path=os.path.abspath(os.path.join(workdir, "outputs")),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=0.1,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            ).distribute(
                protocols=(my_first_protocol, my_second_protocol, my_third_protocol)
            )

            cluster = PyRosettaCluster(
                tasks={
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "ACDEFGHIKLMNPQRSTVWY" * 10,
                },
                input_packed_pose=None,
                seeds=[
                    random.randint(-int((2 ** 32) / 2), int((2 ** 32) / 2) - 1)
                    for _ in range(3)
                ],
                decoy_ids=[random.randint(0, 2) for _ in range(3)],
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=1,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name=None,
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=os.path.abspath(os.path.join(workdir, "outputs")),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=0.5,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )

            cluster.distribute(
                protocols=[my_first_protocol, my_second_protocol, my_third_protocol]
            )


class SaveAllTest(unittest.TestCase):
    def test_save_all(self):
        """Smoke test for PyRosettaCluster usage with the save_all attribute enabled."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )

        _total_tasks = 4
        _total_protocols = 4

        def create_tasks():
            for i in range(_total_tasks):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "LYELL" * (i + 1),
                }

        def my_pyrosetta_protocol(input_packed_pose, **kwargs):
            return input_packed_pose

        with tempfile.TemporaryDirectory() as workdir:
            output_path = os.path.join(workdir, "outputs")
            scorefile_name = "test_save_all.json"
            nstruct = 1
            cluster = PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence("TESTING"),
                seeds=None,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=nstruct,
                dashboard_address=None,
                compressed=True,
                logging_level="CRITICAL",
                scorefile_name=scorefile_name,
                project_name="PyRosettaCluster_SaveAllTest",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=0.5,
                sha1=None,
                dry_run=False,
                save_all=True,
                system_info=None,
                pyrosetta_build=None,
            )
            protocol_args = [my_pyrosetta_protocol] * _total_protocols
            cluster.distribute(*protocol_args)

            with open(os.path.join(output_path, scorefile_name), "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), _total_tasks * _total_protocols)
            _decoy_names = []
            for record in data:
                self.assertDictEqual(record["scores"], {})
                _decoy_name = record["metadata"]["decoy_name"]
                self.assertNotIn(_decoy_name, _decoy_names)
                _decoy_names.append(_decoy_name)
                self.assertEqual(record["instance"]["nstruct"], nstruct)

    def test_save_all_dry_run(self):
        """
        Smoke test for PyRosettaCluster usage with the save_all attribute and
        dry_run attribute enabled.
        """
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        _total_tasks = 2
        _total_protocols = 5

        def create_tasks():
            for i in range(_total_tasks):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "LYELL" * (i + 1),
                }

        def my_pyrosetta_protocol(input_packed_pose, **kwargs):
            return input_packed_pose

        with tempfile.TemporaryDirectory() as workdir:
            output_path = os.path.join(workdir, "outputs")
            scorefile_name = "test_save_all.json"
            decoy_dir_name = "test_decoys"
            logs_dir_name = "test_logs"
            PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence("TESTING"),
                seeds=None,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=2,
                dashboard_address=None,
                compressed=True,
                logging_level="CRITICAL",
                scorefile_name=scorefile_name,
                project_name="PyRosettaCluster_SaveAllTest",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                decoy_dir_name=decoy_dir_name,
                logs_dir_name=logs_dir_name,
                simulation_records_in_scorefile=True,
                ignore_errors=False,
                timeout=0.5,
                sha1=None,
                dry_run=True,
                save_all=True,
                system_info=None,
                pyrosetta_build=None,
            ).distribute(protocols=[my_pyrosetta_protocol] * _total_protocols)

            self.assertFalse(os.path.exists(os.path.join(output_path, scorefile_name)))
            self.assertTrue(os.path.exists(os.path.join(output_path, decoy_dir_name)))
            self.assertTrue(os.path.exists(os.path.join(output_path, logs_dir_name)))
            self.assertEqual(
                os.listdir(os.path.join(output_path, decoy_dir_name)), [],
            )
            self.assertNotEqual(
                os.listdir(os.path.join(output_path, logs_dir_name)), [],
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
