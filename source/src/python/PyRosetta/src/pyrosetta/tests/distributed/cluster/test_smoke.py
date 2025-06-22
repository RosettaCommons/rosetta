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
import numpy
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import random
import subprocess
import sys
import tempfile
import unittest
import warnings

try:
    import cloudpickle
    from dask.distributed import Client, LocalCluster
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_smoke' requires the "
        + "third-party packages 'dask' and 'cloudpickle' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
        + "https://pypi.org/project/cloudpickle/\n"
    )
    raise

from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.utility import get_package_version

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    Serialization,
    get_scores_dict,
    produce,
    reserve_scores,
    run,
    update_scores,
)
from pyrosetta.distributed.cluster.exceptions import WorkerError


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
                max_delay_time=0.5,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )
            cluster = PyRosettaCluster(**instance_kwargs)
            cluster.distribute(
                my_pyrosetta_protocol,
            )

            instance_kwargs.update({"protocols": my_pyrosetta_protocol})
            produce(**instance_kwargs)
            run(**instance_kwargs)

    def test_ignore_errors(self):
        """Test PyRosettaCluster usage with user-provided PyRosetta protocol error."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300",
            set_logging_handler="logging",
        )

        def create_tasks():
            yield {
                "extra_options": "-ex1 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
            }

        def protocol_with_error(packed_pose, **kwargs):
            raise NotImplementedError("Testing an error in a user-provided PyRosetta protocol.")

        with tempfile.TemporaryDirectory() as workdir:
            instance_kwargs = dict(
                tasks=create_tasks,
                input_packed_pose=None,
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
                ignore_errors=True,
                timeout=0.1,
                max_delay_time=3.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )
            cluster = PyRosettaCluster(**instance_kwargs)
            cluster.distribute(protocol_with_error)

        with tempfile.TemporaryDirectory() as workdir:
            instance_kwargs = dict(
                tasks=create_tasks,
                input_packed_pose=None,
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
                max_delay_time=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )
            cluster = PyRosettaCluster(**instance_kwargs)
            with self.assertRaises(WorkerError):
                cluster.distribute(protocol_with_error)


class SmokeTestMulti(unittest.TestCase):
    _ref_kwargs = {
        "test_str": "testing",
        "test_int": 100,
        "test_float": 12345.67890,
        "test_complex": 3j,
        "test_list": list(range(5)),
        "test_tuple": tuple(range(6)),
        "test_range": range(7),
        "test_dict": dict(enumerate(range(8))),
        "test_set": set(range(9)),
        "test_frozenset": frozenset({"foo", "bar"}),
        "test_bool": True,
        "test_bytes": b"Bytes",
        "test_bytearray": bytearray(10),
        "test_none": None,
        "test_memoryview": memoryview(bytes(3)),
    }

    @classmethod
    def setUpClass(cls):
        cloudpickle_version = get_package_version("cloudpickle")
        test_script = os.path.join(os.path.dirname(__file__), "skip_cloudpickle_version.py")
        p = subprocess.run("{0} {1}".format(sys.executable, test_script), shell=True)
        if p.returncode == 0:
            print("Running {0} tests because cloudpickle version {1} can pickle Pose objects.".format(
                    cls.__name__, cloudpickle_version
                )
            )
        elif p.returncode == 1:
            raise unittest.SkipTest(
                "Skipping {0} tests because cloudpickle version {1} cannot pickle Pose objects.".format(
                    cls.__name__, cloudpickle_version
                )
            )
        else:
            raise RuntimeError("Got exit code {0} from running {1}".format(
                    p.returncode, test_script
                )
            )

    def test_smoke_multi(self):
        """Smoke test for PyRosettaCluster usage with multiple protocols."""
        import pyrosetta
        import pyrosetta.distributed
        import pyrosetta.distributed.io as io

        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )

        def _get_xyz_from_packed_pose(packed_pose):
            return list(packed_pose.pose.residue(1).atom(1).xyz())

        def create_tasks():
            for i in range(1, 4):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "TEST" * i,
                    "task_packed_pose": io.pose_from_sequence("PACKED" * i),
                }

        def my_first_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            from pyrosetta.rosetta.core.pose import setPoseExtraScore

            pose = pyrosetta.io.pose_from_sequence(kwargs["seq"])
            setPoseExtraScore(pose, "test_setPoseExtraScore", 123)

            self.assertIn("task_packed_pose", kwargs)
            self.assertIsInstance(kwargs["task_packed_pose"], PackedPose)
            self.assertIn("task_packed_pose", kwargs["PyRosettaCluster_task"])
            self.assertEqual(
                kwargs["PyRosettaCluster_task"]["task_packed_pose"].pose.sequence(),
                kwargs["task_packed_pose"].pose.sequence(),
            )
            saved_packed_pose = io.to_packed(pose.clone())
            kwargs["saved_packed_pose"] = saved_packed_pose
            kwargs["saved_xyz"] = _get_xyz_from_packed_pose(saved_packed_pose)

            saved_pose = pose.clone()
            kwargs["saved_pose"] = saved_pose

            return [pose.clone() for _ in range(2)] + [kwargs]

        @reserve_scores
        def my_second_protocol(packed_pose, **kwargs):
            import pyrosetta  # noqa
            import pyrosetta.distributed.io as io  # noqa

            self.assertEqual(
                dict(packed_pose.scores),
                {**dict(packed_pose.scores), **{"test_setPoseExtraScore": 123}},
            )
            packed_pose.scores.clear()
            self.assertDictEqual({}, packed_pose.scores)
            pose = io.to_pose(packed_pose)

            self.assertIn("task_packed_pose", kwargs)
            task_packed_pose = kwargs["task_packed_pose"]
            self.assertIsInstance(task_packed_pose, PackedPose)

            self.assertIn("saved_packed_pose", kwargs)
            saved_packed_pose = kwargs["saved_packed_pose"]
            self.assertIsInstance(saved_packed_pose, PackedPose)

            self.assertIn("saved_xyz", kwargs)
            saved_xyz = kwargs["saved_xyz"]
            self.assertListEqual(_get_xyz_from_packed_pose(saved_packed_pose), saved_xyz)

            self.assertIn("saved_pose", kwargs)
            saved_pose = kwargs["saved_pose"]
            self.assertIsInstance(saved_pose, Pose)

            if not "saved_variable" in kwargs.keys():
                kwargs["saved_variable"] = {"foo": "bar"}
            else:
                self.assertIn("task_packed_pose", kwargs["PyRosettaCluster_task"])
                kwargs["PyRosettaCluster_task"].pop("task_packed_pose")
                self.assertNotIn(
                    "task_packed_pose",
                    kwargs["PyRosettaCluster_task"],
                    msg="Object of type `PackedPose` is not JSON serializable.",
                )
                self.assertIn(
                    "task_packed_pose",
                    kwargs,
                    msg="`kwargs` with `PackedPose` object(s) should not get JSON serialized.",
                )

            yield kwargs
            yield pose.clone()

        def my_third_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            self.assertEqual(
                dict(packed_pose.scores),
                {**dict(packed_pose.scores), **{"test_setPoseExtraScore": 123}},
            )

            self.assertIn("task_packed_pose", kwargs)
            self.assertIsInstance(kwargs["task_packed_pose"], PackedPose)

            self.assertIn("saved_packed_pose", kwargs)
            saved_packed_pose = kwargs["saved_packed_pose"]
            self.assertIsInstance(saved_packed_pose, PackedPose)

            self.assertIn("saved_xyz", kwargs)
            saved_xyz = kwargs["saved_xyz"]
            self.assertListEqual(_get_xyz_from_packed_pose(saved_packed_pose), saved_xyz)

            self.assertIn("saved_variable", kwargs)
            self.assertDictEqual(kwargs["saved_variable"], {"foo": "bar"})

            return my_second_protocol(packed_pose, **kwargs)

        def my_fourth_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            self.assertFalse(packed_pose.pose.empty())
            self.assertIsInstance(packed_pose, PackedPose)
            self.assertIsInstance(kwargs, dict)
            self.assertNotIn("task_packed_pose", kwargs["PyRosettaCluster_task"])
            self.assertIn("saved_packed_pose", kwargs)
            self.assertIsInstance(kwargs["saved_packed_pose"], PackedPose)

            yield None
            yield kwargs
            yield None

        def my_fifth_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            self.assertIsInstance(packed_pose, PackedPose)
            self.assertTrue(packed_pose.pose.empty())
            self.assertIsInstance(kwargs, dict)
            self.assertIn("saved_packed_pose", kwargs)
            self.assertIsInstance(kwargs["saved_packed_pose"], PackedPose)

            return kwargs

        def my_sixth_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            self.assertIsInstance(packed_pose, PackedPose)
            self.assertTrue(packed_pose.pose.empty())
            self.assertIsInstance(kwargs, dict)
            self.assertNotIn("task_packed_pose", kwargs["PyRosettaCluster_task"])
            self.assertIn("saved_packed_pose", kwargs)
            self.assertIsInstance(kwargs["saved_packed_pose"], PackedPose)

            return None, None, None

        def my_seventh_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            self.assertIsInstance(packed_pose, PackedPose)
            self.assertTrue(packed_pose.pose.empty())
            self.assertIsInstance(kwargs, dict)
            self.assertIn("saved_packed_pose", kwargs)
            self.assertIsInstance(kwargs["saved_packed_pose"], PackedPose)

            return kwargs, None

        def my_eighth_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            kwargs = dict(foo="bar", baz="qux")

            yield kwargs

        def my_ninth_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            self.assertIn("PyRosettaCluster_protocols_container", kwargs)
            self.assertIn("PyRosettaCluster_logging_file", kwargs)
            self.assertIn("PyRosettaCluster_task", kwargs)
            self.assertIn("PyRosettaCluster_protocol_name", kwargs)
            self.assertIn("PyRosettaCluster_protocols", kwargs)
            self.assertIn("PyRosettaCluster_protocol_number", kwargs)
            self.assertIn("PyRosettaCluster_datetime_start", kwargs)
            self.assertIn("PyRosettaCluster_seeds", kwargs)
            self.assertIn("PyRosettaCluster_decoy_ids", kwargs)
            self.assertIn("foo", kwargs)
            self.assertIn("baz", kwargs)

            kwargs = SmokeTestMulti._ref_kwargs

            return kwargs

        def my_tenth_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            ref_kwargs = SmokeTestMulti._ref_kwargs
            for k in ref_kwargs.keys():
                self.assertIn(k, kwargs)
                self.assertEqual(ref_kwargs[k], kwargs[k])

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
                max_delay_time=0.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            ).distribute(
                protocols=(
                    my_first_protocol,
                    my_second_protocol,
                    my_third_protocol,
                    my_fourth_protocol,
                    my_fifth_protocol,
                    my_sixth_protocol,
                    my_seventh_protocol,
                    my_eighth_protocol,
                    my_ninth_protocol,
                    my_tenth_protocol,
                ),
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

            from pyrosetta.rosetta.core.pose import setPoseExtraScore

            pose = pyrosetta.io.pose_from_sequence(kwargs["seq"])
            setPoseExtraScore(pose, "test_setPoseExtraScore", 123)
            self.assertEqual(kwargs["PyRosettaCluster_protocol_number"], 0)
            return [pose.clone() for _ in range(3)]

        @reserve_scores
        def my_second_protocol(packed_pose, **kwargs):
            import pyrosetta  # noqa
            import pyrosetta.distributed.io as io

            self.assertEqual(
                dict(packed_pose.scores),
                {**dict(packed_pose.scores), **{"test_setPoseExtraScore": 123}},
            )
            packed_pose.scores.clear()
            self.assertDictEqual({}, packed_pose.scores)
            self.assertIn(kwargs["PyRosettaCluster_protocol_number"], [1, 2])
            pose = io.to_pose(packed_pose)
            for _ in range(3):
                yield pose.clone()

        def my_third_protocol(packed_pose, **kwargs):
            self.assertEqual(
                dict(packed_pose.scores),
                {**dict(packed_pose.scores), **{"test_setPoseExtraScore": 123}},
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
                seeds=[random.randint(-int((2**32) / 2), int((2**32) / 2) - 1) for _ in range(3)],
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
            ).distribute(protocols=(my_first_protocol, my_second_protocol, my_third_protocol))

            cluster = PyRosettaCluster(
                tasks={
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "ACDEFGHIKLMNPQRSTVWY" * 10,
                },
                input_packed_pose=None,
                seeds=[random.randint(-int((2**32) / 2), int((2**32) / 2) - 1) for _ in range(3)],
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

            cluster.distribute(protocols=[my_first_protocol, my_second_protocol, my_third_protocol])


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
                max_delay_time=1.0,
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
                os.listdir(os.path.join(output_path, decoy_dir_name)),
                [],
            )
            self.assertNotEqual(
                os.listdir(os.path.join(output_path, logs_dir_name)),
                [],
            )


class SerializationTest(unittest.TestCase):
    def test_serialization(self):
        """Smoke test for PyRosettaCluster PackedPose serialization round-trip."""
        for _compression in ("xz", "zlib", "bz2", True, False, None):
            scores = {"test_str": "foo", "test_int": 123, "test_float": numpy.pi}
            for _test_case in range(3):
                input_packed_pose = io.pose_from_sequence("A" * 100)
                if _test_case == 0:
                    # Test scores not cached in Pose
                    input_packed_pose.scores = scores
                elif _test_case == 1:
                    # Test scores not cached in Pose, then update scores
                    input_packed_pose.scores = scores
                    input_packed_pose = update_scores(input_packed_pose)
                elif _test_case == 2:
                    # Test scores cached in Pose
                    input_packed_pose = input_packed_pose.update_scores(scores)

                serializer = Serialization(compression=_compression)
                compressed_packed_pose = serializer.compress_packed_pose(input_packed_pose)
                output_packed_pose = serializer.decompress_packed_pose(compressed_packed_pose)

                _error_msg = f"Failed on test case {_test_case} with compression {_compression}"
                self.assertLess(
                    sys.getsizeof(compressed_packed_pose),
                    sys.getsizeof(input_packed_pose.pickled_pose),
                    msg=_error_msg,
                )
                self.assertLess(
                    sys.getsizeof(compressed_packed_pose),
                    sys.getsizeof(output_packed_pose.pickled_pose),
                    msg=_error_msg,
                )
                if _compression in (False, None):
                    self.assertEqual(id(input_packed_pose), id(output_packed_pose), msg=_error_msg)
                else:
                    self.assertNotEqual(id(input_packed_pose), id(output_packed_pose), msg=_error_msg)
                self.assertEqual(scores, output_packed_pose.scores, msg=_error_msg)
                self.assertSetEqual(
                    set(input_packed_pose.scores.keys()),
                    set(output_packed_pose.scores.keys()),
                    msg=_error_msg,
                )
                for scoretype in input_packed_pose.scores.keys():
                    input_value = input_packed_pose.scores[scoretype]
                    output_value = output_packed_pose.scores[scoretype]
                    if isinstance(input_value, str):
                        self.assertEqual(input_value, output_value, msg=_error_msg)
                    elif isinstance(input_value, (int, float)):
                        self.assertAlmostEqual(input_value, output_value, places=6, msg=_error_msg)

                if _compression not in (False, None):
                    if _test_case == 0:
                        self.assertEqual(scores, input_packed_pose.scores, msg=_error_msg)
                        self.assertEqual(
                            input_packed_pose.scores, output_packed_pose.scores, msg=_error_msg
                        )
                        self.assertNotEqual(
                            input_packed_pose.pickled_pose,
                            output_packed_pose.pickled_pose,
                            msg=_error_msg,
                        )
                    elif _test_case in (1, 2):
                        self.assertEqual(scores, input_packed_pose.scores, msg=_error_msg)
                        self.assertEqual(input_packed_pose.scores, output_packed_pose.scores, msg=_error_msg)
                        self.assertEqual(
                            input_packed_pose.pickled_pose,
                            output_packed_pose.pickled_pose,
                            msg=_error_msg,
                        )


class MultipleClientsTest(unittest.TestCase):
    def test_clients(self):
        """Smoke test for the use case of multiple clients with PyRosettaCluster."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        with tempfile.TemporaryDirectory() as workdir:
            with warnings.catch_warnings():
                # Catch 'ResourceWarning: unclosed <socket.socket ...' from distributed/node.py:235
                # Catch 'UserWarning: Port 8787 is already in use' from distributed/node.py:240
                # Catch 'DeprecationWarning: `np.bool8` is a deprecated alias for `np.bool_`.  (Deprecated NumPy 1.24)' from bokeh/core/property/primitive.py:37
                # Catch 'DeprecationWarning: pkg_resources is deprecated as an API.' from jupyter_server_proxy/config.py:10
                warnings.simplefilter("ignore", category=ResourceWarning)
                warnings.simplefilter("ignore", category=UserWarning)
                warnings.simplefilter("ignore", category=DeprecationWarning)
                cluster_1 = LocalCluster(
                    n_workers=1,
                    threads_per_worker=1,
                    dashboard_address=None,
                    local_directory=workdir,
                )
                cluster_2 = LocalCluster(
                    n_workers=1,
                    threads_per_worker=1,
                    dashboard_address=None,
                    local_directory=workdir,
                )

            client_1 = Client(cluster_1)
            client_1_repr = repr(client_1)

            client_2 = Client(cluster_2)        
            client_2_repr = repr(client_2)

            clients = [client_1, client_2]
            clients_indices = [0, 1, 1, 0]
        
            def create_tasks(client_1_repr=client_1_repr, client_2_repr=client_2_repr):
                for _ in range(1, 7):
                    yield {
                        "extra_options": "-ex1 -multithreading:total_threads 1",
                        "set_logging_handler": "logging",
                        "client_1_repr": client_1_repr,
                        "client_2_repr": client_2_repr,
                    }

            def my_pyrosetta_protocol_1(packed_pose, **kwargs):
                _client_repr = kwargs["PyRosettaCluster_client_repr"]
                _protocol_number = kwargs["PyRosettaCluster_protocol_number"]
                _client_1_repr = kwargs["client_1_repr"]
                _client_2_repr = kwargs["client_2_repr"]

                if _protocol_number == 0:
                    self.assertEqual(_client_repr, _client_1_repr)
                    self.assertNotEqual(_client_repr, _client_2_repr)
                elif _protocol_number == 2:
                    self.assertEqual(_client_repr, _client_2_repr)
                    self.assertNotEqual(_client_repr, _client_1_repr)
                self.assertNotIn(_protocol_number, (1, 3))
            
            def my_pyrosetta_protocol_2(packed_pose, **kwargs):
                _client_repr = kwargs["PyRosettaCluster_client_repr"]
                _protocol_number = kwargs["PyRosettaCluster_protocol_number"]
                _client_1_repr = kwargs["client_1_repr"]
                _client_2_repr = kwargs["client_2_repr"]

                if _protocol_number == 1:
                    self.assertEqual(_client_repr, _client_2_repr)
                    self.assertNotEqual(_client_repr, _client_1_repr)
                elif _protocol_number == 3:
                    self.assertEqual(_client_repr, _client_1_repr)
                    self.assertNotEqual(_client_repr, _client_2_repr)
                self.assertNotIn(_protocol_number, (0, 2))
            
            protocols = [
                my_pyrosetta_protocol_1,
                my_pyrosetta_protocol_2,
                my_pyrosetta_protocol_1,
                my_pyrosetta_protocol_2,
            ]
        
            instance_kwargs = dict(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence("TESTING"),
                seeds=None,
                decoy_ids=None,
                client=None,
                clients=clients,
                clients_indices=clients_indices,
                protocols=protocols,
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
                max_delay_time=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )
            produce(**instance_kwargs)

            for worker in cluster_1.workers.values():
                worker.close_gracefully()
            client_1.close()
            for worker in cluster_2.workers.values():
                worker.close_gracefully()
            client_2.close()


class ResourcesTest(unittest.TestCase):
    def test_resources(self):
        """Smoke test for the use case of abstract resource constraints for dask workers with PyRosettaCluster."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        with tempfile.TemporaryDirectory() as workdir:
            with warnings.catch_warnings():
                # Catch 'ResourceWarning: unclosed <socket.socket ...' from distributed/node.py:235
                # Catch 'UserWarning: Port 8787 is already in use' from distributed/node.py:240
                # Catch 'DeprecationWarning: `np.bool8` is a deprecated alias for `np.bool_`.  (Deprecated NumPy 1.24)' from bokeh/core/property/primitive.py:37
                # Catch 'DeprecationWarning: pkg_resources is deprecated as an API.' from jupyter_server_proxy/config.py:10
                warnings.simplefilter("ignore", category=ResourceWarning)
                warnings.simplefilter("ignore", category=UserWarning)
                warnings.simplefilter("ignore", category=DeprecationWarning)
                cluster = LocalCluster(
                    n_workers=1,
                    threads_per_worker=1,
                    dashboard_address=None,
                    local_directory=workdir,
                    resources={"foo": 1, "bar": 2, "baz": 3},
                )
            client = Client(cluster)
        
            def create_tasks():
                for _ in range(1, 4):
                    yield {
                        "extra_options": "-ex1 -multithreading:total_threads 1",
                        "set_logging_handler": "logging",
                    }

            def my_pyrosetta_protocol_1(packed_pose, **kwargs):   
                return packed_pose.update_scores(sequence=packed_pose.pose.sequence())
            
            def my_pyrosetta_protocol_2(packed_pose, **kwargs):
                return packed_pose
            
            protocols = [
                my_pyrosetta_protocol_1,
                my_pyrosetta_protocol_2,
                my_pyrosetta_protocol_2,
            ]
            resources = [{"foo": 1}, {"bar": 2}, {"baz": 3 - 0.5}]
            output_path = os.path.join(workdir, "outputs")
            decoy_dir_name = "decoys"
            sequence = "TESTING"

            instance_kwargs = dict(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence(sequence),
                seeds=None,
                decoy_ids=None,
                client=client,
                clients=None,
                clients_indices=None,
                protocols=protocols,
                resources=resources,
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
                logging_level="DEBUG",
                scorefile_name=None,
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=False,
                decoy_dir_name=decoy_dir_name,
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=0.1,
                max_delay_time=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )
            produce(**instance_kwargs)

            for worker in cluster.workers.values():
                worker.close_gracefully()
            client.close()

            decoy_files = glob.glob(os.path.join(output_path, decoy_dir_name, "*", "*.bz2"))
            self.assertEqual(len(decoy_files), 3)
            for decoy_file in decoy_files:
                scores_dict = get_scores_dict(decoy_file)
                self.assertEqual(scores_dict["scores"]["sequence"], sequence)

    def test_resources_clients(self):
        """Smoke test for the use case of abstract resource constraints for dask workers with multiple clients in PyRosettaCluster."""
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        with tempfile.TemporaryDirectory() as workdir:
            with warnings.catch_warnings():
                # Catch 'ResourceWarning: unclosed <socket.socket ...' from distributed/node.py:235
                # Catch 'UserWarning: Port 8787 is already in use' from distributed/node.py:240
                # Catch 'DeprecationWarning: `np.bool8` is a deprecated alias for `np.bool_`.  (Deprecated NumPy 1.24)' from bokeh/core/property/primitive.py:37
                # Catch 'DeprecationWarning: pkg_resources is deprecated as an API.' from jupyter_server_proxy/config.py:10
                warnings.simplefilter("ignore", category=ResourceWarning)
                warnings.simplefilter("ignore", category=UserWarning)
                warnings.simplefilter("ignore", category=DeprecationWarning)
                cluster_1 = LocalCluster(
                    n_workers=1,
                    threads_per_worker=1,
                    dashboard_address=None,
                    local_directory=workdir,
                    resources={"FOO": 1},
                )
                cluster_2 = LocalCluster(
                    n_workers=1,
                    threads_per_worker=1,
                    dashboard_address=None,
                    local_directory=workdir,
                    resources={"BAR": 1e9,},
                )

            client_1 = Client(cluster_1)
            client_1_repr = repr(client_1)

            client_2 = Client(cluster_2)        
            client_2_repr = repr(client_2)

            clients = [client_1, client_2]
            clients_indices = [0, 1]
        
            def create_tasks(client_1_repr=client_1_repr, client_2_repr=client_2_repr):
                for _ in range(1, 4):
                    yield {
                        "extra_options": "-ex1 -multithreading:total_threads 1",
                        "set_logging_handler": "logging",
                        "client_1_repr": client_1_repr,
                        "client_2_repr": client_2_repr,
                    }

            def my_pyrosetta_protocol_1(packed_pose, **kwargs):
                _client_repr = kwargs["PyRosettaCluster_client_repr"]
                _protocol_number = kwargs["PyRosettaCluster_protocol_number"]
                self.assertEqual(_protocol_number, 0)
                self.assertNotEqual(_protocol_number, 1)
                _client_1_repr = kwargs["client_1_repr"]
                _client_2_repr = kwargs["client_2_repr"]
                self.assertEqual(_client_repr, _client_1_repr)
                self.assertNotEqual(_client_repr, _client_2_repr)

                return packed_pose.update_scores(sequence=packed_pose.pose.sequence())
                
            def my_pyrosetta_protocol_2(packed_pose, **kwargs):
                _client_repr = kwargs["PyRosettaCluster_client_repr"]
                _protocol_number = kwargs["PyRosettaCluster_protocol_number"]
                self.assertEqual(_protocol_number, 1)
                self.assertNotEqual(_protocol_number, 0)
                _client_1_repr = kwargs["client_1_repr"]
                _client_2_repr = kwargs["client_2_repr"]
                self.assertEqual(_client_repr, _client_2_repr)
                self.assertNotEqual(_client_repr, _client_1_repr)

                return packed_pose
            
            protocols = [my_pyrosetta_protocol_1, my_pyrosetta_protocol_2]
            resources = [{"FOO": 1}, {"BAR": 9e8}]
            output_path = os.path.join(workdir, "outputs")
            decoy_dir_name = "decoys"
            sequence = "TESTING"

            instance_kwargs = dict(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence(sequence),
                seeds=None,
                decoy_ids=None,
                client=None,
                clients=clients,
                clients_indices=clients_indices,
                protocols=protocols,
                resources=resources,
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
                logging_level="DEBUG",
                scorefile_name=None,
                project_name="PyRosettaCluster_Tests",
                simulation_name=None,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=False,
                decoy_dir_name=decoy_dir_name,
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=0.1,
                max_delay_time=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )
            produce(**instance_kwargs)

            decoy_files = glob.glob(os.path.join(output_path, decoy_dir_name, "*", "*.bz2"))
            self.assertEqual(len(decoy_files), 3)
            for decoy_file in decoy_files:
                scores_dict = get_scores_dict(decoy_file)
                self.assertEqual(scores_dict["scores"]["sequence"], sequence)

            for worker in cluster_1.workers.values():
                worker.close_gracefully()
            client_1.close()
            for worker in cluster_2.workers.values():
                worker.close_gracefully()
            client_2.close()

if __name__ == "__main__":
    unittest.main(verbosity=2)
