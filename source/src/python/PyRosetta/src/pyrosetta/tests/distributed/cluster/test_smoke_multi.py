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

import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import random
import subprocess
import sys
import tempfile
import unittest
import uuid

try:
    import cloudpickle # noqa
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_smoke_multi' requires the "
        + "third-party package 'cloudpickle' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/cloudpickle/\n",
        flush=True,
    )
    raise

from pyrosetta import Pose
from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    reserve_scores,
)
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.utility import get_package_version


class SmokeTestMulti(unittest.TestCase):
    """Smoke test for using multiple user-provided PyRosetta protocols in PyRosettaCluster."""

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
        p = subprocess.run(
            "{0} {1}".format(sys.executable, test_script),
            shell=True,
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )
        if p.returncode == 0:
            print("Running {0} tests because cloudpickle version {1} can pickle Pose objects.".format(
                    cls.__name__, cloudpickle_version
                ),
                flush=True,
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

    def tearDown(self):
        sys.stdout.flush()

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
            for i in range(1, 3):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "set_logging_handler": "logging",
                    "seq": "TEST" * i,
                    "task_b64_pose": io.to_base64(io.pose_from_sequence("PACKED" * i)),
                }

        def my_first_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            from pyrosetta.rosetta.core.pose import setPoseExtraScore

            init_options = pyrosetta.get_init_options(compressed=False, as_dict=True)
            self.assertEqual(init_options["multithreading:total_threads"], ["1"])
            self.assertEqual(init_options["packing:ex1"], ["true"])
            self.assertEqual(init_options["packing:ex2aro"], ["true"])

            pose = pyrosetta.io.pose_from_sequence(kwargs["seq"])
            setPoseExtraScore(pose, "test_setPoseExtraScore", 123)

            kwargs["task_packed_pose"] = io.to_packed(kwargs["task_b64_pose"])
            self.assertIn("task_packed_pose", kwargs)
            self.assertIsInstance(kwargs["task_packed_pose"], PackedPose)
            self.assertIn("task_b64_pose", kwargs["PyRosettaCluster_task"])
            self.assertEqual(
                io.to_pose(kwargs["PyRosettaCluster_task"]["task_b64_pose"]).sequence(),
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

            self.assertIsInstance(packed_pose.scores, dict)
            self.assertEqual(packed_pose.scores, {})
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.all_keys)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.extra)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.extra.real)
            self.assertEqual(packed_pose.pose.cache.extra.real["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache.extra["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache["test_setPoseExtraScore"], 123.0)

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
                self.assertIn("task_b64_pose", kwargs["PyRosettaCluster_task"])
                kwargs["PyRosettaCluster_task"].pop("task_b64_pose")
                self.assertNotIn(
                    "task_b64_pose",
                    kwargs["PyRosettaCluster_task"],
                    msg="Object of type `str` is in 'PyRosettaCluster_task' entry.",
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

            self.assertIsInstance(packed_pose.scores, dict)
            self.assertEqual(packed_pose.scores, {})
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.all_keys)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.extra)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.extra.real)
            self.assertEqual(packed_pose.pose.cache.extra.real["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache.extra["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache["test_setPoseExtraScore"], 123.0)

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
            self.assertNotIn("task_b64_pose", kwargs["PyRosettaCluster_task"])
            self.assertIn("saved_packed_pose", kwargs)
            self.assertIsInstance(kwargs["saved_packed_pose"], PackedPose)
            tmp_path = kwargs["PyRosettaCluster_tmp_path"]
            self.assertTrue(os.path.exists(tmp_path))
            protocol_number = kwargs.get("PyRosettaCluster_protocol_number")
            kwargs[f"tmp_path_protocol_{protocol_number}"] = tmp_path

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
            tmp_path = kwargs["PyRosettaCluster_tmp_path"]
            self.assertTrue(os.path.exists(tmp_path))
            protocol_number = kwargs.get("PyRosettaCluster_protocol_number")
            previous_tmp_path_key = f"tmp_path_protocol_{protocol_number - 1}"
            self.assertIn(previous_tmp_path_key, kwargs)
            self.assertNotEqual(kwargs[previous_tmp_path_key], tmp_path)

            return kwargs

        def my_sixth_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            self.assertIsInstance(packed_pose, PackedPose)
            self.assertTrue(packed_pose.pose.empty())
            self.assertIsInstance(kwargs, dict)
            self.assertNotIn("task_b64_pose", kwargs["PyRosettaCluster_task"])
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

            kwargs = dict(foo="bar", baz="qux", PyRosettaCluster_foo="quux")

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
            self.assertIn("PyRosettaCluster_client_repr", kwargs)
            self.assertIn("PyRosettaCluster_output_path", kwargs)
            self.assertIn("PyRosettaCluster_seed", kwargs)
            self.assertIn("PyRosettaCluster_tmp_path", kwargs)
            self.assertIn("foo", kwargs)
            self.assertIn("baz", kwargs)
            self.assertNotIn("PyRosettaCluster_foo", kwargs)

            kwargs = SmokeTestMulti._ref_kwargs

            return kwargs

        def my_tenth_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            ref_kwargs = SmokeTestMulti._ref_kwargs
            for k in ref_kwargs.keys():
                self.assertIn(k, kwargs)
                self.assertEqual(ref_kwargs[k], kwargs[k])

            self.assertNotIn("options", kwargs)
            self.assertNotIn("extra_options", kwargs)
            self.assertNotIn("set_logging_handler", kwargs)
            self.assertNotIn("notebook", kwargs)
            self.assertNotIn("silent", kwargs)
            init_options = pyrosetta.get_init_options(compressed=False, as_dict=True)
            self.assertEqual(init_options["out:levels"], ["all:warning"])
            self.assertEqual(init_options["packing:ex1"], ["true"])
            self.assertEqual(init_options["packing:ex2aro"], ["true"])

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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=os.path.join(workdir, "outputs"),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                filter_results=False,
                norm_task_options=True,
                output_init_file=None,
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

            self.assertIsInstance(packed_pose.scores, dict)
            self.assertEqual(packed_pose.scores, {})
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.all_keys)
            self.assertEqual(packed_pose.pose.cache["test_setPoseExtraScore"], 123.0)
            if kwargs["PyRosettaCluster_protocol_number"] == 1:
                self.assertEqual(kwargs["PyRosettaCluster_protocol_name"], "my_second_protocol")
                self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.extra)
                self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.extra.real)
                self.assertEqual(packed_pose.pose.cache.extra.real["test_setPoseExtraScore"], 123.0)
                self.assertEqual(packed_pose.pose.cache.extra["test_setPoseExtraScore"], 123.0)
            elif kwargs["PyRosettaCluster_protocol_number"] == 2:
                self.assertEqual(kwargs["PyRosettaCluster_protocol_name"], "my_third_protocol")
                self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.metrics)
                self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.metrics.real)
                self.assertEqual(packed_pose.pose.cache.metrics.real["test_setPoseExtraScore"], 123.0)
                self.assertEqual(packed_pose.pose.cache.metrics["test_setPoseExtraScore"], 123.0)
            packed_pose.scores.clear()
            self.assertDictEqual({}, packed_pose.scores)
            self.assertIn(kwargs["PyRosettaCluster_protocol_number"], [1, 2])
            pose = io.to_pose(packed_pose)
            pose.cache.clear() # Clear scoreterms from `pose.cache.extra.real`
            for _ in range(3):
                yield pose.clone()

        def my_third_protocol(packed_pose, **kwargs):
            self.assertIsInstance(packed_pose.scores, dict)
            self.assertEqual(packed_pose.scores, {})
            # `reserve_scores` decorator on 'my_second_protocol' sets scoreterms in `packed_pose.pose.cache.metrics`
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.all_keys)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.metrics)
            self.assertIn("test_setPoseExtraScore", packed_pose.pose.cache.metrics.real)
            self.assertEqual(packed_pose.pose.cache.metrics.real["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache.metrics["test_setPoseExtraScore"], 123.0)
            self.assertEqual(packed_pose.pose.cache["test_setPoseExtraScore"], 123.0)
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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=os.path.abspath(os.path.join(workdir, "outputs")),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                filter_results=True,
                norm_task_options=None,
                output_init_file=None,
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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=os.path.abspath(os.path.join(workdir, "outputs")),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=False,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                filter_results=True,
                norm_task_options=None,
                output_init_file=None,
            )

            cluster.distribute(protocols=[my_first_protocol, my_second_protocol, my_third_protocol])
