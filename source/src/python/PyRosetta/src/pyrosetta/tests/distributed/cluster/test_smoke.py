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


import bz2
import collections
import glob
import inspect
import json
import logging
import numpy
import os
import psutil
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import random
import signal
import subprocess
import sys
import tempfile
import time
import unittest
import uuid
import warnings

try:
    import cloudpickle
    import dask
    import pandas
    import toolz
    from dask.distributed import Client, LocalCluster
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_smoke' requires the "
        + "third-party packages 'dask', 'pandas', and 'cloudpickle' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
        + "https://pypi.org/project/pandas/\n"
        + "https://pypi.org/project/cloudpickle/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

try:
    import cryptography
    has_cryptography = True
except ImportError as ex:
    has_cryptography = False

from pprint import pprint
from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.utility import get_package_version
from pyrosetta.utility.initialization import (
    PyRosettaInitFileReader,
    PyRosettaInitFileSerializer,
    PyRosettaInitFileWriter,
)

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    Serialization,
    export_init_file,
    iterate,
    get_scores_dict,
    produce,
    reserve_scores,
    run,
    update_scores,
)
from pyrosetta.distributed.cluster.exceptions import WorkerError
from pyrosetta.distributed.cluster.init_files import InitFileSigner
from pyrosetta.distributed.cluster.io import (
    METADATA_INPUT_DECOY_KEY,
    METADATA_OUTPUT_DECOY_KEY,
    get_poses_from_init_file,
    secure_read_pickle,
    sign_init_file_metadata_and_poses,
    verify_init_file,
    _is_pandas_object_pyarrow_backed,
)
from pyrosetta.distributed.cluster.multiprocessing import run_protocol, user_protocol
from pyrosetta.distributed.cluster.task_registry import DiskTaskRegistry, MemoryTaskRegistry


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
            security = pyrosetta.distributed.cluster.generate_dask_tls_security(
                os.path.join(workdir, "security_test")
            )
            print(
                "Successfully ran `pyrosetta.distributed.cluster.generate_dask_tls_security()` "
                + f"to generate a dask `Security` object: {security}"
            )
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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=os.path.join(workdir, "outputs"),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                output_decoy_types=[".pdb", ".b64_pose"],
                output_scorefile_types=[".json", ".gz", ".xz"],
                filter_results=True,
                security=security,
                norm_task_options=False,
                output_init_file=None,
            )
            if "pandas" not in pyrosetta.secure_unpickle.get_secure_packages():
                with self.assertRaises(AssertionError):  # output_scorefile_types=[".gz", ...] requires 'pandas' as a secure package
                    PyRosettaCluster(**instance_kwargs)
            pyrosetta.secure_unpickle.add_secure_package("pandas")
            cluster = PyRosettaCluster(**instance_kwargs)
            cluster.distribute(
                my_pyrosetta_protocol,
            )
            instance_kwargs.update({"protocols": my_pyrosetta_protocol, "security": False})
            produce(**instance_kwargs)
            instance_kwargs.update({"security": True})
            if has_cryptography:
                print("Using the installed 'cryptography' package to generate a dask `Security.temporary()` object...")
                run(**instance_kwargs)
            else:
                with self.assertRaises(ImportError) as ex:
                    run(**instance_kwargs)
                print(f"Successfully caught exception: {ex.exception}")

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

        _sep = "*" * 60
        print(f"{_sep} Begin testing PyRosettaCluster(ignore_errors=...) {_sep}")

        with tempfile.TemporaryDirectory() as workdir:
            ignore_errors = True
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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=os.path.join(workdir, "outputs"),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=ignore_errors,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                norm_task_options=True,
                filter_results=False,
                output_init_file=None,
            )
            cluster = PyRosettaCluster(**instance_kwargs)
            cluster.distribute(protocol_with_error)

        with tempfile.TemporaryDirectory() as workdir:
            ignore_errors = False
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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=os.path.join(workdir, "outputs"),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=ignore_errors,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                norm_task_options=True,
                filter_results=False,
                output_init_file=None,
            )
            cluster = PyRosettaCluster(**instance_kwargs)
            with self.assertRaises(WorkerError):
                cluster.distribute(protocol_with_error)

        print(f"{_sep} End testing PyRosettaCluster(ignore_errors=...) {_sep}")


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
            for i in range(1, 3):
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
            self.assertNotIn("task_packed_pose", kwargs["PyRosettaCluster_task"])
            self.assertIn("saved_packed_pose", kwargs)
            self.assertIsInstance(kwargs["saved_packed_pose"], PackedPose)
            tmp_path = kwargs["PyRosettaCluster_tmp_path"]
            self.assertTrue(os.path.exists(tmp_path))

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
            author = "Username"
            email = "test@example"
            license = "LICENSE.PyRosetta.md"
            project_name = None
            simulation_name = "SaveAllTest"
            input_packed_pose = io.pose_from_sequence("TESTING")
            output_decoy_types = [".pdb", ".init"]
            output_scorefile_types = [".json"]
            compressed = True
            cluster = PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=input_packed_pose,
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
                compressed=compressed,
                logging_level="CRITICAL",
                scorefile_name=scorefile_name,
                project_name=project_name,
                simulation_name=simulation_name,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=True,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                sha1=None,
                dry_run=False,
                save_all=True,
                system_info=None,
                pyrosetta_build=None,
                author=author,
                email=email,
                license=license,
                filter_results=True,
                output_decoy_types=output_decoy_types,
                output_scorefile_types=output_scorefile_types,
                norm_task_options=None,
            ) # Test 'output_init_file' default
            protocol_args = [my_pyrosetta_protocol] * _total_protocols
            cluster.distribute(*protocol_args)

            with open(os.path.join(output_path, scorefile_name), "r") as f:
                data = [json.loads(line) for line in f]
            self.assertEqual(len(data), _total_tasks * _total_protocols)

            _project_name = "PyRosettaCluster" if project_name is None else project_name
            _simulation_name = "PyRosettaCluster" if simulation_name is None else simulation_name
            init_file = os.path.join(output_path, f"{_project_name}_{_simulation_name}_pyrosetta.init")
            if compressed:
                init_file += ".bz2"
            self.assertTrue(os.path.isfile(init_file))
            for record in data:
                self.assertTrue(os.path.samefile(init_file, record["metadata"]["init_file"]))
            if compressed:
                init_data = io.read_init_file(init_file)
            else:
                init_data = PyRosettaInitFileReader.read_json(init_file)

            # Verify output PyRosetta initialization file
            self.assertEqual(init_data["author"], author)
            self.assertEqual(init_data["email"], email)
            self.assertEqual(init_data["license"], license)
            self.assertEqual(init_data["pyrosetta_build"], pyrosetta._build_signature())
            self.assertListEqual(init_data["options"]["run:constant_seed"], ["true"])

            # Verify metadata in output PyRosetta initialization file
            metadata = init_data["metadata"]
            self.assertIsInstance(metadata, dict)
            self.assertEqual(metadata["comment"], "Generated by PyRosettaCluster")
            self.assertEqual(metadata["version"], pyrosetta.distributed.cluster.__version__)
            if isinstance(input_packed_pose, (Pose, PackedPose)):
                self.assertIn(METADATA_INPUT_DECOY_KEY, metadata)
            else:
                self.assertNotIn(METADATA_INPUT_DECOY_KEY, metadata)
            self.assertNotIn(METADATA_OUTPUT_DECOY_KEY, metadata)

            self.assertIn("sha256", metadata)
            sha256 = metadata.pop("sha256", None)
            self.assertIn("signature", metadata)
            signature = metadata.pop("signature", None)
            signer = InitFileSigner(
                input_packed_pose=input_packed_pose,
                output_packed_pose=None,
                metadata=metadata,
            )
            self.assertTrue(signer.verify_sha256(sha256))
            self.assertTrue(signer.verify_signature(signature))

            # Export PyRosetta initialization file with output decoy
            record = data[0] # Use first output decoy
            output_file = record["metadata"]["output_file"]

            if ".init" in output_decoy_types:
                if compressed:
                    decoy_init_file = os.path.splitext(os.path.splitext(output_file)[0])[0] + ".init.bz2"
                    self.assertTrue(os.path.isfile(decoy_init_file))
                    with open(decoy_init_file, "rb") as fbz2:
                        decoy_init_data = PyRosettaInitFileReader.from_json(bz2.decompress(fbz2.read()).decode())
                    _decoy_init_data = io.read_init_file(decoy_init_file)
                    self.assertDictEqual(decoy_init_data, _decoy_init_data)
                else:
                    decoy_init_file = os.path.splitext(output_file)[0] + ".init"
                    self.assertTrue(os.path.isfile(decoy_init_file))
                    decoy_init_data = PyRosettaInitFileReader.read_json(decoy_init_file)
                    _decoy_init_data = io.read_init_file(decoy_init_file)
                    self.assertDictEqual(decoy_init_data, _decoy_init_data)

                # Verify decoy output PyRosetta initialization file
                self.assertEqual(decoy_init_data["author"], author)
                self.assertEqual(decoy_init_data["email"], email)
                self.assertEqual(decoy_init_data["license"], license)
                self.assertEqual(decoy_init_data["pyrosetta_build"], pyrosetta._build_signature())
                self.assertListEqual(decoy_init_data["options"]["run:constant_seed"], ["true"])

                # Verify metadata in decoy output PyRosetta initialization file
                decoy_metadata = decoy_init_data["metadata"]
                self.assertIsInstance(decoy_metadata, dict)
                self.assertEqual(decoy_metadata["comment"], "Generated by PyRosettaCluster")
                self.assertEqual(decoy_metadata["version"], pyrosetta.distributed.cluster.__version__)
                if isinstance(input_packed_pose, (Pose, PackedPose)):
                    self.assertIn(METADATA_INPUT_DECOY_KEY, decoy_metadata)
                else:
                    self.assertNotIn(METADATA_INPUT_DECOY_KEY, decoy_metadata)
                self.assertIn(METADATA_OUTPUT_DECOY_KEY, decoy_metadata)
                poses_output_idx = decoy_metadata[METADATA_OUTPUT_DECOY_KEY]
                self.assertEqual(poses_output_idx, 0)

                self.assertIn("sha256", decoy_metadata)
                decoy_sha256 = decoy_metadata.pop("sha256")
                self.assertIn("signature", decoy_metadata)
                decoy_signature = decoy_metadata.pop("signature")
                decoy_signer_1 = InitFileSigner(
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=io.to_packed(io.to_pose(decoy_init_data["poses"][poses_output_idx])),
                    metadata=decoy_metadata,
                )
                self.assertTrue(decoy_signer_1.verify_sha256(decoy_sha256))
                self.assertTrue(decoy_signer_1.verify_signature(decoy_signature))
                decoy_signer_2 = InitFileSigner(
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=io.pose_from_init_file(decoy_init_file),
                    metadata=decoy_metadata,
                )
                self.assertTrue(decoy_signer_2.verify_sha256(decoy_sha256))
                self.assertTrue(decoy_signer_2.verify_signature(decoy_signature))
                decoy_signer_3 = InitFileSigner(
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=io.poses_from_init_file(decoy_init_file)[poses_output_idx],
                    metadata=decoy_metadata,
                )
                self.assertTrue(decoy_signer_3.verify_sha256(decoy_sha256))
                self.assertTrue(decoy_signer_3.verify_signature(decoy_signature))

            custom_output_init_file = os.path.join(workdir, "custom_output.init")
            export_init_file(
                output_file,
                output_init_file=custom_output_init_file,  # Test custom output file path
                compressed=False,  # Test not compressed
            )
            self.assertTrue(os.path.isfile(custom_output_init_file))
            export_init_file(
                output_file,
                output_init_file=custom_output_init_file,  # Test custom output file path
                compressed=True,  # Test compressed
            )
            custom_output_init_bz2_file = custom_output_init_file + ".bz2"
            self.assertTrue(os.path.isfile(custom_output_init_bz2_file))
            if ".init" not in output_decoy_types:
                export_init_file(
                    output_file,
                    output_init_file=None,  # Test default output file path
                    compressed=compressed,
                )
                exported_init_file = (
                    output_file[: -len(f"{output_decoy_types[0]}.bz2")]
                    if compressed
                    else output_file[: -len(output_decoy_types[0])]
                ) + (".init.bz2" if compressed else ".init")
            else:
                exported_init_file = (
                    output_file[: -len(f"{output_decoy_types[0]}.bz2")]
                    if compressed
                    else output_file[: -len(output_decoy_types[0])]
                ) + "_my_custom.init"
                export_init_file(
                    output_file,
                    output_init_file=exported_init_file,
                    compressed=compressed,
                )
                if compressed:
                    exported_init_file += ".bz2"
            self.assertTrue(os.path.isfile(exported_init_file))
            exported_init_data = io.read_init_file(exported_init_file)

            # Verify exported PyRosetta initialization file
            self.assertEqual(exported_init_data["author"], author)
            self.assertEqual(exported_init_data["email"], email)
            self.assertEqual(exported_init_data["license"], license)
            self.assertEqual(exported_init_data["pyrosetta_build"], pyrosetta._build_signature())
            self.assertListEqual(exported_init_data["options"]["run:constant_seed"], ["true"])

            # Verify metadata in exported PyRosetta initialization file
            exported_metadata = exported_init_data["metadata"]
            self.assertIsInstance(exported_metadata, dict)
            self.assertEqual(exported_metadata["comment"], "Generated by PyRosettaCluster")
            self.assertEqual(exported_metadata["version"], pyrosetta.distributed.cluster.__version__)
            if isinstance(input_packed_pose, (Pose, PackedPose)):
                self.assertIn(METADATA_INPUT_DECOY_KEY, exported_metadata)
            else:
                self.assertNotIn(METADATA_INPUT_DECOY_KEY, exported_metadata)
            self.assertIn(METADATA_OUTPUT_DECOY_KEY, exported_metadata)
            self.assertIn("sha256", exported_metadata)
            exported_sha256 = exported_metadata.pop("sha256", None)
            self.assertIn("signature", exported_metadata)
            exported_signature = exported_metadata.pop("signature", None)
            exported_signer = InitFileSigner(
                input_packed_pose=input_packed_pose,
                output_packed_pose=io.pose_from_init_file(exported_init_file), # Output decoy is first object in 'poses' key value
                metadata=exported_metadata,
            )
            self.assertTrue(exported_signer.verify_sha256(exported_sha256))
            self.assertTrue(exported_signer.verify_signature(exported_signature))

            _decoy_names = set()
            for record in data:
                self.assertEqual(init_data["pyrosetta_build"], record["instance"]["pyrosetta_build"])
                self.assertEqual(exported_init_data["pyrosetta_build"], record["instance"]["pyrosetta_build"])
                for key in ("author", "email", "license"):
                    self.assertIn(key, record["instance"])
                    self.assertNotIn(key, record["metadata"])
                    self.assertNotIn(key, record["scores"])
                self.assertEqual(record["instance"]["author"], author)
                self.assertEqual(record["instance"]["email"], email)
                self.assertEqual(record["instance"]["license"], license)
                self.assertNotIn("init_file", record["instance"])
                self.assertIn("init_file", record["metadata"])
                self.assertTrue(os.path.samefile(record["metadata"]["init_file"], init_file))
                self.assertDictEqual(record["scores"], {})
                _decoy_name = record["metadata"]["decoy_name"]
                self.assertNotIn(_decoy_name, _decoy_names)
                _decoy_names.add(_decoy_name)
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
        _nstruct = 2

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
            init_file = os.path.join(output_path, "test.init")
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
                nstruct=_nstruct,
                dashboard_address=None,
                compressed=True,
                logging_level="CRITICAL",
                scorefile_name=scorefile_name,
                project_name="PyRosettaCluster_SaveAllTest",
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=output_path,
                decoy_dir_name=decoy_dir_name,
                logs_dir_name=logs_dir_name,
                simulation_records_in_scorefile=True,
                ignore_errors=True,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=True,
                save_all=True,
                system_info=None,
                pyrosetta_build=None,
                filter_results=True,
                norm_task_options=None,
                output_init_file=init_file,
            ).distribute(protocols=[my_pyrosetta_protocol] * _total_protocols)

            self.assertFalse(os.path.exists(os.path.join(output_path, scorefile_name)))
            self.assertFalse(os.path.exists(init_file))
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

                if _compression in (False, None):
                    # Test PackedPose size
                    if _test_case == 0:
                        self.assertGreater(
                            sys.getsizeof(compressed_packed_pose.pickled_pose),
                            sys.getsizeof(input_packed_pose.pickled_pose),
                            msg=_error_msg,
                        ) # `PackedPose.scores` get cached in `Pose.cache`
                    elif _test_case in (1, 2):
                        self.assertEqual(
                            sys.getsizeof(compressed_packed_pose.pickled_pose),
                            sys.getsizeof(input_packed_pose.pickled_pose),
                            msg=_error_msg,
                        ) # `PackedPose.scores` are already cached in `Pose.cache`
                    self.assertEqual(
                        sys.getsizeof(compressed_packed_pose.pickled_pose),
                        sys.getsizeof(output_packed_pose.pickled_pose),
                        msg=_error_msg,
                    ) # `PackedPose.scores` are already cached in `Pose.cache`
                    # Test PackedPose identity
                    self.assertEqual(id(compressed_packed_pose), id(output_packed_pose), msg=_error_msg)
                else:
                    # Test PackedPose size
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
                    # Test PackedPose identity
                    self.assertNotEqual(id(compressed_packed_pose), id(output_packed_pose), msg=_error_msg)
                # Test PackedPose identity
                self.assertNotEqual(id(input_packed_pose), id(compressed_packed_pose), msg=_error_msg)
                self.assertNotEqual(id(input_packed_pose), id(output_packed_pose), msg=_error_msg)
                # Test PackedPose scores
                if _test_case == 0: # Compression has no effect when user doesn't cache scores in Pose object
                    self.assertNotEqual(input_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.scores, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.pose.cache, {}, msg=_error_msg)
                    self.assertEqual(output_packed_pose.scores, {}, msg=_error_msg)
                    self.assertEqual(output_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.scores, output_packed_pose.pose.cache, msg=_error_msg)
                    self.assertEqual(output_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertSetEqual(
                        set(input_packed_pose.scores.keys()),
                        set(output_packed_pose.pose.cache.all_keys),
                        msg=_error_msg,
                    )
                    for scoretype in input_packed_pose.scores.keys():
                        input_value = input_packed_pose.scores[scoretype]
                        output_value = output_packed_pose.pose.cache[scoretype]
                        self.assertEqual(input_value, output_value, msg=_error_msg)
                    self.assertNotEqual(
                        input_packed_pose.pickled_pose,
                        output_packed_pose.pickled_pose,
                        msg=_error_msg,
                    )
                elif _test_case in (1, 2): # Compression has no effect when user caches scores in Pose object (preferred syntax)
                    self.assertNotEqual(input_packed_pose.scores, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.scores, {}, msg=_error_msg)
                    self.assertEqual(input_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertEqual(output_packed_pose.scores, {}, msg=_error_msg)
                    self.assertEqual(input_packed_pose.scores, output_packed_pose.scores, msg=_error_msg)
                    self.assertEqual(output_packed_pose.pose.cache, scores, msg=_error_msg)
                    self.assertEqual(input_packed_pose.pose.cache, output_packed_pose.pose.cache, msg=_error_msg)
                    self.assertSetEqual(
                        set(input_packed_pose.pose.cache.all_keys),
                        set(output_packed_pose.pose.cache.all_keys),
                        msg=_error_msg,
                    )
                    for scoretype in input_packed_pose.pose.cache.all_keys:
                        input_value = input_packed_pose.pose.cache[scoretype]
                        output_value = output_packed_pose.pose.cache[scoretype]
                        self.assertEqual(input_value, output_value, msg=_error_msg)
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
                for _ in range(1, 3):
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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=os.path.join(workdir, "outputs"),
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                filter_results=False,
                norm_task_options=None,
                output_init_file=None,
            )
            produce(**instance_kwargs)

            client_1.close()
            cluster_1.close()
            client_2.close()
            cluster_2.close()


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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=False,
                decoy_dir_name=decoy_dir_name,
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                filter_results=True,
                norm_task_options=None,
                output_init_file=None,
            )
            produce(**instance_kwargs)

            client.close()
            cluster.close()

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
                simulation_name=uuid.uuid4().hex,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=False,
                decoy_dir_name=decoy_dir_name,
                logs_dir_name="logs",
                ignore_errors=True,
                timeout=1.0,
                max_delay_time=3.0,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
                filter_results=True,
                norm_task_options=None,
                output_init_file=None,
            )
            produce(**instance_kwargs)

            decoy_files = glob.glob(os.path.join(output_path, decoy_dir_name, "*", "*.bz2"))
            self.assertEqual(len(decoy_files), 3)
            for decoy_file in decoy_files:
                scores_dict = get_scores_dict(decoy_file)
                self.assertEqual(scores_dict["scores"]["sequence"], sequence)

            client_1.close()
            cluster_1.close()
            client_2.close()
            cluster_2.close()


class ScoresTest(unittest.TestCase):
    _value = 1e1

    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        cls.input_packed_pose = io.pose_from_sequence("TEST")
        cls.workdir = tempfile.TemporaryDirectory()
        cls.decoy_dir_name = "decoys"
        cls.instance_kwargs = dict(
            tasks=ScoresTest.create_task,
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
            compressed=True,
            logging_level="INFO",
            scorefile_name=None,
            project_name="PyRosettaCluster_Tests",
            simulation_name=uuid.uuid4().hex,
            environment=None,
            simulation_records_in_scorefile=False,
            decoy_dir_name=cls.decoy_dir_name,
            logs_dir_name="logs",
            ignore_errors=True,
            timeout=1.0,
            max_delay_time=3.0,
            sha1=None,
            dry_run=False,
            save_all=False,
            system_info=None,
            pyrosetta_build=None,
            filter_results=True,
            norm_task_options=None,
            output_init_file=None,
        )

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    @staticmethod
    def create_task():
        yield {
            "extra_options": "-ex1 -multithreading:total_threads 1",
            "set_logging_handler": "logging",
        }

    @staticmethod
    def identity_protocol(packed_pose, **kwargs):
        import pyrosetta
        import pyrosetta.distributed.io as io

        return packed_pose

    @staticmethod
    @reserve_scores
    def reserved_scores_protocol(packed_pose, **kwargs):
        import pyrosetta
        import pyrosetta.distributed.io as io

        pose = packed_pose.pose
        pose.cache.clear()
        packed_pose = io.to_packed(pose)
        packed_pose.scores.clear()

        return packed_pose

    @staticmethod
    def add_detached_scores_protocol(packed_pose, **kwargs):
        import pyrosetta
        import pyrosetta.distributed.io as io

        packed_pose = packed_pose.update_scores(attached_score=ScoresTest._value)
        packed_pose.scores["detached_score"] = ScoresTest._value

        return packed_pose

    @staticmethod
    def protocol_with_secure_package_pandas(packed_pose, **kwargs):
        assert "pandas" not in pyrosetta.secure_unpickle.get_secure_packages()
        pyrosetta.secure_unpickle.add_secure_package("pandas")
        assert "pandas" in pyrosetta.secure_unpickle.get_secure_packages()
        if "secure_packages" in kwargs and "pyarrow" in kwargs["secure_packages"]:
            assert "pyarrow" not in pyrosetta.secure_unpickle.get_secure_packages()
            pyrosetta.secure_unpickle.add_secure_package("pyarrow")
            assert "pyarrow" in pyrosetta.secure_unpickle.get_secure_packages()
        _ = packed_pose.pose.cache["df"]
        return packed_pose

    @staticmethod
    def protocol_without_secure_package_pandas(packed_pose, **kwargs):
        _ = packed_pose.pose.cache["df"]
        return packed_pose

    def get_scores_dict(self, output_path):
        decoy_files = glob.glob(os.path.join(output_path, self.decoy_dir_name, "*", "*.bz2"))
        self.assertEqual(len(decoy_files), 1)
        scores_dict = get_scores_dict(next(iter(decoy_files)))

        return scores_dict

    def setup_input_packed_pose(self):
        pose = io.to_pose(self.input_packed_pose).clone()
        pose.cache.clear()
        input_packed_pose = io.to_packed(pose)
        input_packed_pose.scores.clear()

        return input_packed_pose

    def test_detached_scores(self):
        """Test saving detached scores in PyRosettaCluster with/without compression."""
        for compression in (True, False):
            input_packed_pose = self.setup_input_packed_pose()
            input_packed_pose = self.input_packed_pose.update_scores(attached_score=ScoresTest._value)
            input_packed_pose.scores["detached_score"] = ScoresTest._value
            output_path = os.path.join(self.workdir.name, f"test_detached_scores_{compression}")
            run(
                **{
                    **self.instance_kwargs,
                    "input_packed_pose": input_packed_pose,
                    "protocols": ScoresTest.identity_protocol,
                    "compression": compression,
                    "output_path": output_path,
                }
            )
            scores_dict = self.get_scores_dict(output_path)
            for key in ("attached_score", "detached_score"):
                self.assertIn(
                    key,
                    scores_dict["scores"],
                    msg=f"Saving score '{key}' failed with compression={compression}",
                )
                self.assertEqual(scores_dict["scores"][key], ScoresTest._value)

    def test_detached_scores_with_reserve_scores(self):
        """Test saving detached scores in PyRosettaCluster with/without compression with `reserve_scores` decorator."""
        for compression in (True, False):
            input_packed_pose = self.setup_input_packed_pose()
            input_packed_pose = self.input_packed_pose.update_scores(attached_score=ScoresTest._value)
            input_packed_pose.scores["detached_score"] = ScoresTest._value
            output_path = os.path.join(self.workdir.name, f"test_detached_scores_with_reserve_scores_{compression}")
            run(
                **{
                    **self.instance_kwargs,
                    "input_packed_pose": input_packed_pose,
                    "protocols": ScoresTest.reserved_scores_protocol,
                    "compression": compression,
                    "output_path": output_path,
                }
            )
            scores_dict = self.get_scores_dict(output_path)
            for key in ("attached_score", "detached_score"):
                self.assertIn(
                    key,
                    scores_dict["scores"],
                    msg=f"Saving score '{key}' failed with compression={compression}",
                )
                self.assertEqual(scores_dict["scores"][key], ScoresTest._value)

    def test_detached_scores_in_protocol(self):
        """Test saving detached scores in PyRosettaCluster protocol with/without compression."""
        for compression in (True, False):
            input_packed_pose = self.setup_input_packed_pose()
            output_path = os.path.join(self.workdir.name, f"test_detached_scores_in_protocol_{compression}")
            run(
                **{
                    **self.instance_kwargs,
                    "input_packed_pose": input_packed_pose,
                    "protocols": ScoresTest.add_detached_scores_protocol,
                    "compression": compression,
                    "output_path": output_path,
                }
            )
            scores_dict = self.get_scores_dict(output_path)
            for key in ("attached_score", "detached_score"):
                self.assertIn(
                    key,
                    scores_dict["scores"],
                    msg=f"Saving score '{key}' failed with compression={compression}",
                )
                self.assertEqual(scores_dict["scores"][key], ScoresTest._value)

    def test_secure_packages_billiard(self):
        """
        Test caching a `pandas.DataFrame` with and without adding 'pandas'
        as a secure package in the billiard subprocess.
        """
        pyrosetta.secure_unpickle.add_secure_package("pandas")
        df = pandas.DataFrame().from_dict({0: ["foo"], 1: ["bar"]})
        if _is_pandas_object_pyarrow_backed(df):
            # If the cached `pandas.DataFrame` object uses Arrow-backed dtypes, then
            # PyRosetta requires 'pyarrow' to be in the unpickle-allowed list during
            # output decoy parsing. Certain `pandas` versions (with Python-3.13+)
            # use Arrow-backed dtypes for `pandas.DataFrame` objects by default.
            pyrosetta.secure_unpickle.add_secure_package("pyarrow")
        input_pose = self.input_packed_pose.pose.clone()
        input_pose.cache["df"] = df # Cache `pandas.DataFrame` object
        secure_packages = pyrosetta.secure_unpickle.get_secure_packages()
        # Test a protocol that does not add 'pandas' to the unpickle-allowed list,
        # and does not access the cached `pandas.DataFrame`; this tests that
        # PyRosettaCluster infrastructure does not trigger deserialization alone
        run(
            **{
                **self.instance_kwargs,
                "tasks": [{**task, "secure_packages": secure_packages} for task in ScoresTest.create_task()],
                "input_packed_pose": input_pose.clone(),
                "output_path": os.path.join(self.workdir.name, "test_secure_packages_billiard_1"),
                "ignore_errors": False,
                "protocols": ScoresTest.identity_protocol,
            }
        )
        # Test a protocol that adds 'pandas' to the unpickle-allowed list, and
        # accesses the cached `pandas.DataFrame`; this tests that the billiard
        # subprocess requires adding 'pandas' to the unpickle-allowed list
        # before data access, even though the client process has already added it
        run(
            **{
                **self.instance_kwargs,
                "tasks": [{**task, "secure_packages": secure_packages} for task in ScoresTest.create_task()],
                "input_packed_pose": input_pose.clone(),
                "output_path": os.path.join(self.workdir.name, "test_secure_packages_billiard_2"),
                "ignore_errors": False,
                "protocols": ScoresTest.protocol_with_secure_package_pandas,
            }
        )
        _sep = "*" * 60
        print(f"{_sep} Begin testing expected `UnpickleSecurityError` in billiard subprocess {_sep}")
        with self.assertRaises(WorkerError):
            # Test a protocol that does not add 'pandas' to the unpickle-allowed list,
            # and then accesses the cached `pandas.DataFrame`; this tests that the
            # billiard subprocess requires adding 'pandas' to the unpickle-allowed list
            # before data access, even though the client process has already added it,
            # leading to an intentionally raised `WorkerError` exception
            run(
                **{
                    **self.instance_kwargs,
                    "tasks": [{**task, "secure_packages": secure_packages} for task in ScoresTest.create_task()],
                    "input_packed_pose": input_pose.clone(),
                    "output_path": os.path.join(self.workdir.name, "test_secure_packages_billiard_3"),
                    "ignore_errors": False,
                    "protocols": ScoresTest.protocol_without_secure_package_pandas,
                }
            )
        print(f"{_sep} End testing expected `UnpickleSecurityError` in billiard subprocess {_sep}")


class TestInitFileSigner(unittest.TestCase):
    def test_init_file_signer(self):
        with tempfile.TemporaryDirectory() as workdir:
            output_init_file = os.path.join(workdir, "pyrosetta.init")
            input_packed_pose = io.pose_from_sequence("START")
            output_packed_pose = io.pose_from_sequence("FINAL")
            output_pose = output_packed_pose.pose
            pyrosetta.rosetta.core.pose.add_comment(
                output_pose,
                "REMARK PyRosettaCluster:",
                "{}",
            )
            output_packed_pose = io.to_packed(output_pose)
            poses = [output_packed_pose, input_packed_pose, input_packed_pose.clone()]
            metadata = {"comment": "foo", METADATA_INPUT_DECOY_KEY: "bar", METADATA_OUTPUT_DECOY_KEY: 10, "sha256": "test", "signature": 123}
            with self.assertRaises(ValueError):  # Fails 'sha256' and 'signature' verification
                verify_init_file(
                    output_init_file,
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=output_packed_pose,
                    metadata=metadata,
                )

            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                verbose=False,
            )
            with self.assertRaises(TypeError):  # Fails because metadata `METADATA_INPUT_DECOY_KEY` value is a `str` object
                _ = get_poses_from_init_file(output_init_file, verify=False)
            _poses = io.poses_from_init_file(output_init_file)
            self.assertEqual(len(_poses), 3)
            self.assertEqual(_poses[0].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertEqual(_poses[1].pose.sequence(), input_packed_pose.pose.sequence())

            metadata = {"comment": "bar", METADATA_INPUT_DECOY_KEY: -5, METADATA_OUTPUT_DECOY_KEY: 10, "sha256": None, "signature": None}
            verify_init_file(
                output_init_file,
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
                metadata=metadata,
            )  # Bypass verification with 'sha256' and 'signature' as `NoneType`
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            with self.assertRaises(IndexError):  # Fails since metadata `METADATA_INPUT_DECOY_KEY` and `METADATA_OUTPUT_DECOY_KEY` values are out of range
                _ = get_poses_from_init_file(output_init_file, verify=False)

            metadata = {"comment": "baz", METADATA_INPUT_DECOY_KEY: 1, METADATA_OUTPUT_DECOY_KEY: 1, "sha256": 256}
            with self.assertRaises(ValueError):  # Fails SHA256 verification
                verify_init_file(
                    output_init_file,
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=output_packed_pose,
                    metadata=metadata,
                )
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            _poses = get_poses_from_init_file(output_init_file, verify=False)  # Skip verification test
            self.assertEqual(len(_poses), 2)
            self.assertEqual(_poses[0].pose.sequence(), input_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[0].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertEqual(_poses[1].pose.sequence(), input_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[1].pose.sequence(), output_packed_pose.pose.sequence())

            metadata = {"comment": "baz", METADATA_INPUT_DECOY_KEY: 0, METADATA_OUTPUT_DECOY_KEY: 0}
            verify_init_file(
                output_init_file,
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
                metadata=metadata,
            )  # Bypass verification with 'sha256' and 'signature' missing
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            _poses = get_poses_from_init_file(output_init_file, verify=False)
            self.assertEqual(len(_poses), 2)
            self.assertEqual(_poses[0].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[0].pose.sequence(), input_packed_pose.pose.sequence())
            self.assertEqual(_poses[1].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[1].pose.sequence(), input_packed_pose.pose.sequence())

            metadata = {"comment": "baz", METADATA_INPUT_DECOY_KEY: 1, METADATA_OUTPUT_DECOY_KEY: 0, "signature": "test"}
            with self.assertRaises(ValueError):  # Fails signature verification
                verify_init_file(
                    output_init_file,
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=output_packed_pose,
                    metadata=metadata,
                )
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            with self.assertRaises(ValueError):
                _ = get_poses_from_init_file(output_init_file, verify=True)  # Fails verification
            _poses = get_poses_from_init_file(output_init_file, verify=False)  # Skip verification
            self.assertEqual(len(_poses), 2)

            self.assertEqual(_poses[0].pose.sequence(), input_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[0].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertEqual(_poses[1].pose.sequence(), output_packed_pose.pose.sequence())
            self.assertNotEqual(_poses[1].pose.sequence(), input_packed_pose.pose.sequence())

            metadata = {
                "comment": "Generated by PyRosettaCluster",
                METADATA_INPUT_DECOY_KEY: 1,
                METADATA_OUTPUT_DECOY_KEY: 0,
                "version": pyrosetta.distributed.cluster.__version__,
            }
            verify_init_file(
                output_init_file,
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
                metadata=metadata,
            )  # Passes signature verification
            pyrosetta.dump_init_file(
                output_init_file,
                poses=poses.copy(),
                metadata=metadata,
                overwrite=True,
                verbose=False,
            )
            _poses = get_poses_from_init_file(output_init_file, verify=True)  # Passes verification

            metadata_signed, poses_signed = sign_init_file_metadata_and_poses(
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
            )
            init_dict = PyRosettaInitFileReader.read_json(output_init_file)
            init_dict["metadata"] = metadata_signed
            init_dict["poses"] = [io.to_base64(p) for p in poses_signed]
            init_dict["md5"] = PyRosettaInitFileSerializer.get_md5(init_dict)
            signed_init_file = os.path.join(workdir, "custom_pyrosetta.init")
            PyRosettaInitFileWriter.write_json(init_dict, signed_init_file)
            verify_init_file(
                output_init_file,
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
                metadata=dict(
                    filter(lambda kv: kv[0] not in ("sha256", "signature"), metadata_signed.items())
                ),
            )  # Passes verification

            for _metadata in [
                "string",
                b"Bytes",
                numpy.pi,
                dict(enumerate(range(20))),
                metadata,
                metadata_signed,
            ]:
                if isinstance(_metadata, bytes):
                    with self.assertRaises(TypeError):  # Object of type bytes is not JSON serializable
                        signer = InitFileSigner(
                            input_packed_pose=input_packed_pose,
                            output_packed_pose=output_packed_pose,
                            metadata=_metadata,
                        )
                else:
                    signer = InitFileSigner(
                        input_packed_pose=input_packed_pose,
                        output_packed_pose=output_packed_pose,
                        metadata=_metadata,
                    )
                signatures_dict = signer.sign()
                self.assertIsInstance(signatures_dict, dict)
                signer.verify_sha256(signatures_dict.get("sha256"))
                signer.verify_signature(signatures_dict.get("signature"))
                signer.verify(*signatures_dict.values())


class TestBase:
    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        cls.input_packed_pose = io.pose_from_sequence("TESTING")
        cls.local_directory = tempfile.TemporaryDirectory()
        cls.local_directory_1 = tempfile.TemporaryDirectory()
        cls.local_directory_2 = tempfile.TemporaryDirectory()
        with warnings.catch_warnings():
            # Catch 'ResourceWarning: unclosed <socket.socket ...' from distributed/node.py:235
            # Catch 'UserWarning: Port 8787 is already in use' from distributed/node.py:240
            # Catch 'DeprecationWarning: `np.bool8` is a deprecated alias for `np.bool_`.  (Deprecated NumPy 1.24)' from bokeh/core/property/primitive.py:37
            # Catch 'DeprecationWarning: pkg_resources is deprecated as an API.' from jupyter_server_proxy/config.py:10
            warnings.simplefilter("ignore", category=ResourceWarning)
            warnings.simplefilter("ignore", category=UserWarning)
            warnings.simplefilter("ignore", category=DeprecationWarning)
            default_cluster = LocalCluster(
                n_workers=1,
                threads_per_worker=1,
                dashboard_address=None,
                local_directory=cls.local_directory.name,
                resources={"CPU": 1},
            )
            cluster_1 = LocalCluster(
                n_workers=1,
                threads_per_worker=1,
                dashboard_address=None,
                local_directory=cls.local_directory_1.name,
                resources={"FOO": 1, "BAZ": 2},
            )
            cluster_2 = LocalCluster(
                n_workers=1,
                threads_per_worker=1,
                dashboard_address=None,
                local_directory=cls.local_directory_2.name,
                resources={"BAR": 1e9, "BAZ": 2},
            )
        cls.default_client = Client(default_cluster)
        cls.clients = [Client(cluster_1), Client(cluster_2)]

    @classmethod
    def tearDownClass(cls):
        cls.local_directory.cleanup()
        cls.local_directory_1.cleanup()
        cls.local_directory_2.cleanup()
        for client in [cls.default_client, *cls.clients]:
            for worker in client.cluster.workers.values():
                worker.close_gracefully()
            client.shutdown()
            client.cluster.close()

    def setUp(self):
        self.workdir = tempfile.TemporaryDirectory()
        self.instance_kwargs = dict(
            seeds=None,
            decoy_ids=None,
            client=None,
            clients=None,
            scheduler=None,
            scratch_dir=os.path.join(self.workdir.name, "scratch"),
            cores=None,
            processes=None,
            memory=None,
            min_workers=1,
            max_workers=1,
            nstruct=1,
            dashboard_address=None,
            compression=True,
            compressed=True,
            logging_level="INFO",
            scorefile_name=None,
            project_name="PyRosettaCluster_Tests",
            simulation_name=uuid.uuid4().hex,
            environment=None,
            output_path=os.path.join(self.workdir.name, "outputs_default"),
            simulation_records_in_scorefile=False,
            decoy_dir_name="test_decoy_dir",
            logs_dir_name="logs",
            ignore_errors=True,
            timeout=1.0,
            max_delay_time=3.0,
            sha1=None,
            system_info=None,
            pyrosetta_build=None,
            dry_run=False,
            save_all=False,
            filter_results=True,
            norm_task_options=None,
            output_init_file=None,
        )

    def tearDown(self):
        self.workdir.cleanup()


class GeneratorTest(TestBase, unittest.TestCase):
    """Smoke test for the use case of the `PyRosettaCluster().generate()` method."""
    _n_tasks = 2
    _n_output_packed_poses = 2
    _parameters = (0.0, 100.0) # `float` objects for `packed_pose.update_scores` values
    _pyrosetta_kwargs = {
        "options": "-mute all",
        "extra_options": "-ex1 -multithreading:total_threads 1",
        "set_logging_handler": "logging",
        }

    @staticmethod
    def parameter_to_str(parameter):
        return f"stored_{int(parameter)}"

    @staticmethod
    def create_tasks(parameter):
        for _ in range(GeneratorTest._n_tasks):
            yield {
                "extra_options": "-ex1 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
                "n_output_packed_poses": GeneratorTest._n_output_packed_poses,
                "parameter": parameter,
            }

    @staticmethod
    def my_pyrosetta_protocol_1(packed_pose, **kwargs):
        value = 1.0
        kwargs["kwargs_key_1"] = value
        packed_pose = packed_pose.update_scores(scores_key_1=value)
        yield packed_pose
        yield kwargs

    @staticmethod
    def my_pyrosetta_protocol_2(packed_pose, **kwargs):
        assert "kwargs_key_1" in kwargs
        assert "scores_key_1" in list(packed_pose.pose.scores)
        value = 2.0
        kwargs["kwargs_key_2"] = value
        parameter = kwargs.get("parameter", "string")
        packed_pose = packed_pose.update_scores(
            {
                "scores_key_2": value,
                "parameter": parameter,
                GeneratorTest.parameter_to_str(parameter): "test",
            }
        )
        packed_poses = [packed_pose.pose.clone() for _ in range(kwargs["n_output_packed_poses"])]
        return (*packed_poses, kwargs)

    @classmethod
    def get_protocols(cls):
        return [cls.my_pyrosetta_protocol_1, cls.my_pyrosetta_protocol_2]

    def test_generate_builtin_clients(self):
        """Test for `PyRosettaCluster().generate()` using built-in client instantiations."""
        instance_kwargs_update = {
            **self.instance_kwargs,
            "tasks": self.create_tasks(parameter=GeneratorTest._parameters[0]),
            "input_packed_pose": self.input_packed_pose,
            "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_builtin_clients_1_{uuid.uuid4().hex}"),
        }
        protocols = self.get_protocols()
        clients_indices = instance_kwargs_update.pop("clients_indices", None)
        resources = instance_kwargs_update.pop("resources", None)
        instance = PyRosettaCluster(**instance_kwargs_update)

        output_packed_pose_results = []
        output_kwargs_results = []
        for output_packed_pose, output_kwargs in instance.generate(
            protocols=protocols,
            clients_indices=clients_indices,
            resources=resources,
        ):
            self.assertIsInstance(output_packed_pose, PackedPose)
            self.assertIsInstance(output_kwargs, dict)
            self.assertIn("kwargs_key_1", output_kwargs)
            self.assertIn("kwargs_key_2", output_kwargs)
            self.assertIn("scores_key_1", list(output_packed_pose.pose.scores))
            self.assertIn("scores_key_2", list(output_packed_pose.pose.scores))
            instance_kwargs_update = {
                **self.instance_kwargs,
                "tasks": self.create_tasks(parameter=GeneratorTest._parameters[1]),
                "input_packed_pose": output_packed_pose,
                "protocols": protocols,
                "clients_indices": clients_indices,
                "resources": resources,
                "simulation_name": "test_generate_builtin_clients",
                "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_builtin_clients_2_{uuid.uuid4().hex}"),
            }
            for output_packed_pose, output_kwargs in iterate(**instance_kwargs_update):
                self.assertIsInstance(output_packed_pose, PackedPose)
                self.assertIsInstance(output_kwargs, dict)
                self.assertIn("kwargs_key_1", output_kwargs)
                self.assertIn("kwargs_key_2", output_kwargs)
                self.assertIn("scores_key_1", list(output_packed_pose.pose.scores))
                self.assertIn("scores_key_2", list(output_packed_pose.pose.scores))
                output_packed_pose_results.append(output_packed_pose)
                output_kwargs_results.append(output_kwargs)
        _n_results_per_parameter = (GeneratorTest._n_tasks * GeneratorTest._n_output_packed_poses)
        _n_output_results = _n_results_per_parameter ** 2
        self.assertEqual(len(output_packed_pose_results), _n_output_results)
        self.assertEqual(len(output_kwargs_results), _n_output_results)
        for packed_pose in output_packed_pose_results:
            self.assertEqual(
                packed_pose.pose.scores.get("parameter", None),
                GeneratorTest._parameters[1],
                msg="Output packed pose does not have the correct value for the 'parameter' scores key."
            )
        for kwargs in output_kwargs_results:
            self.assertEqual(
                kwargs["PyRosettaCluster_task"].get("parameter", None),
                GeneratorTest._parameters[1],
                msg="Output kwargs do not have the correct value for the 'parameter' task key."
            )
        for _parameter in GeneratorTest._parameters:
            scores_key = self.parameter_to_str(_parameter)
            self.assertEqual(
                len(list(filter(lambda packed_pose: scores_key in packed_pose.pose.scores, output_packed_pose_results))),
                _n_output_results,
                msg=f"Packed poses do not have the correct values for the '{scores_key}' scores key."
            )

    def test_generate_user_client(self):
        """
        Test for `PyRosettaCluster().generate()` with `save_all=True`,
        `dry_run=True`, and using a user-provided client.
        """
        instance_kwargs = {
            **self.instance_kwargs,
            "tasks": self.create_tasks(parameter=GeneratorTest._parameters[0]),
            "input_packed_pose": self.input_packed_pose,
            "protocols": self.get_protocols(),
            "client": self.default_client,
            "clients": None,
            "clients_indices": None,
            "resources": None,
            "save_all": True,
            "dry_run": True,
            "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_user_client_1_{uuid.uuid4().hex}"),
        }
        results = []
        for output_packed_pose, _ in iterate(**instance_kwargs):
            instance_kwargs_update = {
                **instance_kwargs,
                "input_packed_pose": output_packed_pose,
                "tasks": self.create_tasks(parameter=GeneratorTest._parameters[1]),
                "client": self.default_client, # Test passing in same client
                "simulation_name": "test_generate_user_client",
                "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_user_client_2_{uuid.uuid4().hex}"),
            }
            for result in iterate(**instance_kwargs_update):
                results.append(result)
        _n_output_packed_poses_save_all = (GeneratorTest._n_output_packed_poses + 1) # Plus one from my_pyrosetta_protocol_1
        _n_results_per_parameter = (GeneratorTest._n_tasks * _n_output_packed_poses_save_all)
        _n_output_results = _n_results_per_parameter ** 2
        self.assertEqual(
            len(results), _n_output_results, msg="Number of results with save_all failed."
        )
        self.assertListEqual(
            os.listdir(os.path.join(instance_kwargs["output_path"], instance_kwargs["decoy_dir_name"])),
            [],
            msg="Dry run failed while yielding results.",
        )
        self.assertFalse(
            os.path.isfile(os.path.join(instance_kwargs["output_path"], "scores.json")),
            msg="Dry run failed while yielding results.",
        )

    def test_generate_multi_user_clients(self):
        """
        Test for `PyRosettaCluster().generate()` with `save_all=True`,
        `dry_run=True`, and using multiple user-provided clients.
        """
        clients_indices = [0, 1]
        resources = [{"FOO": 1}, {"BAR": 9e8}]
        instance_kwargs = {
            **self.instance_kwargs,
            "tasks": self.create_tasks(parameter=GeneratorTest._parameters[0]),
            "input_packed_pose": self.input_packed_pose,
            "protocols": self.get_protocols(),
            "client": None,
            "clients": self.clients,
            "clients_indices": clients_indices,
            "resources": resources,
            "save_all": True,
            "dry_run": True,
            "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_multi_user_clients_1_{uuid.uuid4().hex}"),
        }
        results = []
        for output_packed_pose, _ in iterate(**instance_kwargs):
            instance_kwargs_update = {
                **instance_kwargs,
                "input_packed_pose": output_packed_pose,
                "tasks": self.create_tasks(parameter=GeneratorTest._parameters[1]),
                "client": None,
                "clients": self.clients, # Test passing in same clients
                "clients_indices": clients_indices,
                "resources": resources,
                "simulation_name": "test_generate_multi_user_clients",
                "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_multi_user_clients_2_{uuid.uuid4().hex}"),
            }
            for result in iterate(**instance_kwargs_update):
                results.append(result)
        _n_output_packed_poses_save_all = (GeneratorTest._n_output_packed_poses + 1) # Plus one from my_pyrosetta_protocol_1
        _n_results_per_parameter = (GeneratorTest._n_tasks * _n_output_packed_poses_save_all)
        _n_output_results = _n_results_per_parameter ** 2
        self.assertEqual(
            len(results), _n_output_results, msg="Number of results with save_all failed."
        )
        self.assertListEqual(
            os.listdir(os.path.join(instance_kwargs["output_path"], instance_kwargs["decoy_dir_name"])),
            [],
            msg="Dry run failed while yielding results.",
        )
        self.assertFalse(
            os.path.isfile(os.path.join(instance_kwargs["output_path"], "scores.json")),
            msg="Dry run failed while yielding results.",
        )

    def test_generate_partition_clients(self):
        """
        Test for `PyRosettaCluster().generate()` with `save_all=True`,
        `dry_run=True`, and partitioning user-provided clients over iterations.
        """
        resources = [{"BAZ": 1}, {"BAZ": 2}]
        instance_kwargs = {
            **self.instance_kwargs,
            "tasks": self.create_tasks(parameter=GeneratorTest._parameters[0]),
            "input_packed_pose": self.input_packed_pose,
            "protocols": self.get_protocols(),
            "client": self.clients[0], # Test passing in first client
            "clients": None,
            "clients_indices": None,
            "resources": resources,
            "save_all": True,
            "dry_run": True,
            "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_partition_clients_1_{uuid.uuid4().hex}"),
        }
        results = []
        for output_packed_pose, _ in iterate(**instance_kwargs):
            instance_kwargs_update = {
                **instance_kwargs,
                "input_packed_pose": output_packed_pose,
                "tasks": self.create_tasks(parameter=GeneratorTest._parameters[1]),
                "client": self.clients[1], # Test passing in second client
                "clients": None,
                "clients_indices": None,
                "resources": resources,
                "simulation_name": "test_generate_partition_clients",
                "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_partition_clients_2_{uuid.uuid4().hex}"),
            }
            for result in iterate(**instance_kwargs_update):
                results.append(result)
        _n_output_packed_poses_save_all = (GeneratorTest._n_output_packed_poses + 1) # Plus one from my_pyrosetta_protocol_1
        _n_results_per_parameter = (GeneratorTest._n_tasks * _n_output_packed_poses_save_all)
        _n_output_results = _n_results_per_parameter ** 2
        self.assertEqual(
            len(results), _n_output_results, msg="Number of results with save_all failed."
        )
        self.assertListEqual(
            os.listdir(os.path.join(instance_kwargs["output_path"], instance_kwargs["decoy_dir_name"])),
            [],
            msg="Dry run failed while yielding results.",
        )
        self.assertFalse(
            os.path.isfile(os.path.join(instance_kwargs["output_path"], "scores.json")),
            msg="Dry run failed while yielding results.",
        )


class RuntimeTestLoggingFilter(logging.Filter):
    matches = (
        "with_lock",
        "dry_run",
        "save_all",
        "Percent Complete",
        "simulation complete",
        "Attempted to determine the residue type set of an empty pose",
    )
    def filter(self, record):
        msg = record.getMessage()
        return all(map(lambda s: s not in msg, self.matches))


@unittest.skip("Auxiliary tests for runtime testing.")
class RuntimeTest(TestBase, unittest.TestCase):
    @staticmethod
    def create_simple_tasks(n_tasks=10):
        for i in range(n_tasks):
            yield {
                **GeneratorTest._pyrosetta_kwargs,
                "task": i,
            }

    @staticmethod
    def timing_protocol_1(packed_pose, **kwargs):
        import logging
        _logger = logging.getLogger(kwargs['PyRosettaCluster_protocol_name'])
        _logger.info(f"Running procotol {kwargs['PyRosettaCluster_protocol_name']} with task {kwargs['task']}")
        yield packed_pose

    @staticmethod
    def timing_protocol_2(packed_pose, **kwargs):
        return RuntimeTest.timing_protocol_1(packed_pose, **kwargs)

    @classmethod
    def setup_logger(cls, stream=True):
        _logger = logging.getLogger(__name__)
        if stream:
            _stream_handler = logging.StreamHandler(sys.stdout)
            _logger.addHandler(_stream_handler)
        else:
            _stream_handler = None

        _root_logger = logging.getLogger("root")
        _filter = RuntimeTestLoggingFilter()
        _root_logger.addFilter(_filter)

        _rosetta_logger = logging.getLogger("rosetta")
        _rosetta_logger.addFilter(_filter)

        _distributed_logger = logging.getLogger("pyrosetta.distributed")
        _distributed_logger.addFilter(_filter)

        return _logger, _stream_handler

    @classmethod
    def tear_down_logger(cls, _logger, _stream_handler):
        if _stream_handler is not None:
            _logger.removeHandler(_stream_handler)

    @staticmethod
    def get_mean_dt(ts):
        return numpy.mean([ts[i + 1] - ts[i] for i in range(len(ts) - 1)])

    @staticmethod
    def get_dt(t, ts):
        return t if len(ts) == 0 else (t - ts[-1])

    @unittest.skip("Auxiliary test for runtime testing.")
    def test_timing_single_instance(self):
        """Runtime test with a single PyRosettaCluster instance."""
        _logger, _stream_handler = self.setup_logger(stream=True)
        _logger.info(f"Starting single PyRosettaCluster instance runtime test.")
        prc = PyRosettaCluster(
            **{
                **self.instance_kwargs,
                "input_packed_pose": self.input_packed_pose,
                "tasks": self.create_simple_tasks(),
                "clients": self.clients,
                "dry_run": True,
                "save_all": True,
                "max_delay_time": 0.0,
            }
        )
        prc_iterable = prc.generate(
            self.timing_protocol_1,
            self.timing_protocol_2,
            clients_indices=[0, 1],
            resources=[{"FOO": 1}, {"BAZ": 2}],
        )
        ts = []
        t0 = time.time()
        for i, (_, output_kwargs) in enumerate(prc_iterable):
            t = time.time() - t0
            dt = self.get_dt(t, ts)
            ts.append(t)
            task = output_kwargs.get("task")
            client_repr = output_kwargs.get("PyRosettaCluster_client_repr")
            _logger.info(f"Finished iteration {(i,)} with task {task} on client {client_repr} in {dt:0.3f} seconds")
        mean_dt = self.get_mean_dt(ts)
        _logger.info(f"Average iteration time with a single PyRosettaCluster instance: {mean_dt:0.7f} seconds")
        self.tear_down_logger(_logger, _stream_handler)

    @unittest.skip("Auxiliary test for runtime testing.")
    def test_timing_multi_instance(self):
        """Runtime test for two PyRosettaCluster instances asynchronously generating results."""
        _logger, _stream_handler = self.setup_logger(stream=True)
        _logger.info(f"Starting multiple PyRosettaCluster instance runtime test.")
        setup_kwargs = {
            "dry_run": True,
            "save_all": False,
            "max_delay_time": 0.0,
        }
        instance_kwargs_1 = {
            **self.instance_kwargs,
            **setup_kwargs,
            "input_packed_pose": self.input_packed_pose,
            "tasks": self.create_simple_tasks(),
            "client": self.clients[0],
            "protocols": self.timing_protocol_1,
            "resources": [{"BAZ": 1}],
        }
        prc_iterate = iterate(**instance_kwargs_1)
        instance_kwargs_2 = {
            **self.instance_kwargs,
            **setup_kwargs,
            "client": self.clients[1],
            "simulation_name": "test_timing_multi_instance",
        }
        for k in ("protocols", "clients_indices", "resources"):
            instance_kwargs_2.pop(k, None)
        prc = PyRosettaCluster(**instance_kwargs_2)
        ts = []
        t0 = time.time()
        for i, (output_packed_pose, output_kwargs) in enumerate(prc_iterate):
            t = time.time() - t0
            dt = self.get_dt(t, ts)
            ts.append(t)
            task = output_kwargs.get("task")
            client_repr = output_kwargs.get("PyRosettaCluster_client_repr")
            _logger.info(f"Finished iteration {(i,)} with task {task} on client {client_repr} in {dt:0.3f} seconds")
            prc.tasks = [{**GeneratorTest._pyrosetta_kwargs, "task": task}]
            prc.input_packed_pose = output_packed_pose
            for j, (_, output_kwargs) in enumerate(prc.generate(self.timing_protocol_2, resources=[{"BAR": 1e8}])):
                t = time.time() - t0
                dt = self.get_dt(t, ts)
                ts.append(t)
                task = output_kwargs.get("task")
                client_repr = output_kwargs.get("PyRosettaCluster_client_repr")
                _logger.info(f"Finished iteration {(i, j)} with task {task} on client {client_repr} in {dt:0.3f} seconds")
        mean_dt = self.get_mean_dt(ts)
        _logger.info(f"Average iteration time with multiple PyRosettaCluster instances: {mean_dt:0.7f} seconds")
        self.tear_down_logger(_logger, _stream_handler)


class PrioritiesTest(unittest.TestCase):
    def test_priorities(self):
        """Smoke test for the use case of priorities of tasks with PyRosettaCluster."""
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
                    n_workers=1, # One worker for all tasks to test priorities
                    threads_per_worker=1, # One thread per worker to block recurrently running tasks
                    dashboard_address=None,
                    local_directory=workdir,
                )

            client_1 = Client(cluster_1)

            _n_tasks: int = 10
            _n_protocols: int = 3

            def create_tasks():
                t0 = time.time()
                for i in range(_n_tasks):
                    yield {
                        "extra_options": "-ex1 -multithreading:total_threads 1",
                        "set_logging_handler": "logging",
                        "task": i,
                        "t0": t0,
                    }

            def my_pyrosetta_protocol(packed_pose, **kwargs):
                protocol_number = kwargs["PyRosettaCluster_protocol_number"]
                task = kwargs["task"]
                t0 = kwargs["t0"]
                t1 = time.time()
                dt = round(t1 - t0, 1)
                scoretype = f"task_{task}_protocol_{protocol_number}"
                return packed_pose.update_scores({scoretype: dt})

            protocols = [my_pyrosetta_protocol] * _n_protocols
            clients_indices = [0] * _n_protocols
            scorefile_name = "test_priorities.json"
            decoy_dir_name = "decoys"
            sequence = "TASK/CHAIN"
            # Depth-first task chains (increases priority per recursion, so task chains run to completion)
            depth_first_priorities = [
                list(range(0, _n_protocols * 10, 10)),
                list(range(-_n_protocols * 2, 0, 2)),
            ]
            # Breadth-first task chains (first-in, first-out)
            breadth_first_priorities = [
                None,
                [0] * _n_protocols,
                [25] * _n_protocols,
            ]
            priorities_test_cases = {}
            depth_first_test_cases = []
            breadth_first_test_cases = []
            test_case = 0
            # Register depth-first tests explicitly
            for priorities in depth_first_priorities:
                priorities_test_cases[test_case] = priorities
                depth_first_test_cases.append(test_case)
                test_case += 1
            # Register breadth-first tests explicitly
            for priorities in breadth_first_priorities:
                priorities_test_cases[test_case] = priorities
                breadth_first_test_cases.append(test_case)
                test_case += 1
            # Run test cases
            for test_case, priorities in priorities_test_cases.items():
                if isinstance(priorities, (tuple, list)):
                    output_path = os.path.join(workdir, "outputs" + "_".join(map(str, priorities)))
                else:
                    output_path = os.path.join(workdir, "outputs" + f"_{priorities}")
                instance_kwargs = dict(
                    tasks=create_tasks,
                    input_packed_pose=io.pose_from_sequence(sequence),
                    seeds=None,
                    decoy_ids=None,
                    client=None,
                    clients=[client_1],
                    clients_indices=clients_indices,
                    protocols=protocols,
                    priorities=priorities,
                    resources=None,
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
                    scorefile_name=scorefile_name,
                    project_name="PyRosettaCluster_Tests",
                    simulation_name=uuid.uuid4().hex,
                    environment=None,
                    output_path=output_path,
                    simulation_records_in_scorefile=True,
                    decoy_dir_name=decoy_dir_name,
                    logs_dir_name="logs",
                    ignore_errors=True,
                    timeout=0.5,
                    max_delay_time=0.5,
                    sha1=None,
                    dry_run=False,
                    save_all=False,
                    system_info=None,
                    pyrosetta_build=None,
                    filter_results=True,
                    norm_task_options=None,
                    output_init_file=None,
                    output_decoy_types=[".b64_pose"],
                    output_scorefile_types=[".json"],
                    security=False,
                    max_nonce=99999,
                )
                produce(**instance_kwargs)

                scorefile_path = os.path.join(output_path, scorefile_name)
                # Analyze adjacent protocol transition rate
                results = {}
                with open(scorefile_path, "r") as f:
                    for line in f:
                        d = json.loads(line)
                        results.update(d.get("scores"))
                # Sort result
                sorted_results = sorted(results.items(), key=lambda kv: kv[1])
                # Extract task numbers
                task_positions = collections.defaultdict(list)
                for i, (name, _) in enumerate(sorted_results):
                    task_num = int(name.split("_")[1])
                    task_positions[task_num].append(i)
                # Compute variance per task
                task_variances = [numpy.var(pos) for pos in task_positions.values()]
                # Compute mean variance
                spread_metric = float(numpy.mean(task_variances))
                # Compute theoretical maximum variance
                max_var = (len(sorted_results) - 1) ** 2 / 4
                # Compute normalized variance metric
                norm_variance = spread_metric / max_var
                # Empirically, depth-first runs yield norm variance ~0.01-0.08 and breadth-first ~0.16-0.34.
                # Therefore, the midpoint (~0.12) separates the two behaviors
                norm_variance_threshold = 0.12 # Empirically determined
                _msg = ", ".join(
                    [
                        f"num_protocols={_n_protocols}",
                        f"num_tasks={_n_tasks}",
                        f"priorities={priorities}",
                        f"task_positions={dict(task_positions)}",
                        f"spread_metric={spread_metric}",
                        f"max_var={max_var}",
                        f"norm_variance={norm_variance}",
                        f"norm_variance_threshold={norm_variance_threshold}",
                    ]
                )
                if test_case in depth_first_test_cases:
                    self.assertLess(norm_variance, norm_variance_threshold, msg=_msg)
                elif test_case in breadth_first_test_cases:
                    self.assertGreater(norm_variance, norm_variance_threshold, msg=_msg)
                else:
                    raise NotImplementedError(test_case)
                print(f"Successfully tested `priorities` keyword argument: {_msg}")

            client_1.close()
            cluster_1.close()


class IntentionalError(RuntimeError):
    def __init__(self, *args):
        super().__init__(*args)


class RetriesTest(unittest.TestCase):
    """Smoke tests for the use case of task retries with PyRosettaCluster."""
    _n_tasks: int = 3
    _n_protocols: int = 2
    _n_retries_per_protocol: int = 2  # 3 attempts total
    _count_marker: str = "."
    _counter_filename_template: str = "protocol_number-{protocol_number}__task-{task}__retries.txt"
    _sep = "*" * 60

    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        cls.workdir = tempfile.TemporaryDirectory()
        with warnings.catch_warnings():
            # Catch 'ResourceWarning: unclosed <socket.socket ...' from distributed/node.py:235
            # Catch 'UserWarning: Port 8787 is already in use' from distributed/node.py:240
            # Catch 'DeprecationWarning: `np.bool8` is a deprecated alias for `np.bool_`.  (Deprecated NumPy 1.24)' from bokeh/core/property/primitive.py:37
            # Catch 'DeprecationWarning: pkg_resources is deprecated as an API.' from jupyter_server_proxy/config.py:10
            warnings.simplefilter("ignore", category=ResourceWarning)
            warnings.simplefilter("ignore", category=UserWarning)
            warnings.simplefilter("ignore", category=DeprecationWarning)
            cls.cluster_1 = LocalCluster(
                n_workers=1, # Must be 1 for `RetriesTest` unit tests, otherwise tasks run concurrently and test logic will fail
                threads_per_worker=1, # Must be 1 for `RetriesTest` unit tests, otherwise tasks run concurrently and test logic will fail
                dashboard_address=None,
                local_directory=cls.workdir.name,
            )
        cls.client_1 = Client(cls.cluster_1)
        cls.default_instance_kwargs = dict(
            tasks=RetriesTest.create_tasks,
            input_packed_pose=io.pose_from_sequence("TEST/TASK/RETRIES"),
            seeds=None,
            decoy_ids=None,
            client=None,
            clients=[cls.client_1],
            scheduler=None,
            scratch_dir=cls.workdir.name,
            cores=None,
            processes=None,
            memory=None,
            min_workers=1,
            max_workers=1,
            nstruct=1,
            dashboard_address=None,
            compressed=True,
            logging_level="WARNING",
            scorefile_name="scores.json",
            project_name="PyRosettaCluster_Tests",
            simulation_name=uuid.uuid4().hex,
            environment=None,
            simulation_records_in_scorefile=True,
            decoy_dir_name="decoys",
            logs_dir_name="logs",
            ignore_errors=False, # Never ignore errors for `RetriesTest` test case
            timeout=0.5,
            max_delay_time=0.5,
            sha1=None,
            dry_run=False,
            save_all=False,
            system_info=None,
            pyrosetta_build=None,
            filter_results=True, # Protocols return non-empty PackedPose when succeeding
            norm_task_options=None,
            output_init_file=None,
            output_decoy_types=[".pdb"],
            output_scorefile_types=[".json"],
            security=False,
            max_nonce=99999,
        )
        print(f"{RetriesTest._sep} Begin testing PyRosettaCluster().distribute(retries=...) {RetriesTest._sep}")

    @classmethod
    def tearDownClass(cls):
        # Note: During dask client/cluster shutdown, dask worker processes may still have in-flight tasks
        # or scheduled retries, which can lead to emitted warnings like:
        #   - "distributed.worker.state_machine - WARNING - Async instruction for <Task cancelled ...> ended with CancelledError"
        #   - "UserWarning: semaphore_tracker: There appear to be ... leaked semaphores to clean up at shutdown"
        # These warnings are expected from the `RetriesTest` test case, and do not indicate unit test failures.
        cls.client_1.close()
        cls.cluster_1.close()
        time.sleep(3) # Allow logging messages from worker processes to flush
        cls.workdir.cleanup()
        print(f"{RetriesTest._sep} End testing PyRosettaCluster().distribute(retries=...) {RetriesTest._sep}")

    @staticmethod
    def create_tasks():
        for i in range(RetriesTest._n_tasks):
            yield {
                "extra_options": "-ex1 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
                "task": i,
                "max_retries": RetriesTest._n_retries_per_protocol,
                "count_marker": RetriesTest._count_marker,
                "counter_filename_template": RetriesTest._counter_filename_template,
            }

    @staticmethod
    def attempt_counter_protocol(**kwargs):
        output_path = kwargs["PyRosettaCluster_output_path"]
        protocol_number = kwargs["PyRosettaCluster_protocol_number"]
        task = kwargs["task"]
        max_retries = kwargs["max_retries"]
        max_attempts = max_retries + 1 # Plus one, since retries only occur after the initial attempt
        count_marker = kwargs["count_marker"]
        counter_filename_template = kwargs["counter_filename_template"]
        counter_file = os.path.join(
            output_path,
            counter_filename_template.format(protocol_number=protocol_number, task=task),
        )
        # Log new attempt
        with open(counter_file, "a") as f:
            total_attempts = f.write(count_marker)
        # Get total attempts, including this attempt
        with open(counter_file, "r") as f:
            total_attempts = f.read().count(count_marker)

        return task, total_attempts, max_attempts

    @staticmethod
    def my_protocol(packed_pose, **kwargs):
        return packed_pose, kwargs

    @staticmethod
    def protocol_with_intentional_error(packed_pose, **kwargs):
        task, total_attempts, max_attempts = RetriesTest.attempt_counter_protocol(**kwargs)
        # Always raise error
        raise IntentionalError(f"Task {task} intentionally raised error on attempt {total_attempts}/{max_attempts}")

    @staticmethod
    def protocol_succeeds_on_last_retry(packed_pose, **kwargs):
        task, total_attempts, max_attempts = RetriesTest.attempt_counter_protocol(**kwargs)
        if total_attempts == max_attempts: # Succeed on last retry
            print(f"Protocol `protocol_succeeds_on_last_retry` succeeding on last attempt ({total_attempts}/{max_attempts})!")
            return packed_pose
        else: # Fail if not on last retry
            raise IntentionalError(f"Task {task} intentionally raised error on attempt {total_attempts}/{max_attempts}")

    def test_retries_persistent_errors(self):
        """Test retries of protocols with persistent errors."""
        protocols = [RetriesTest.protocol_with_intentional_error] * RetriesTest._n_protocols
        output_path = os.path.join(self.workdir.name, "outputs_persistent_errors")
        prc = PyRosettaCluster(
            **self.default_instance_kwargs,
            **{"output_path": output_path},
        )
        with self.assertRaises(WorkerError):
            prc.distribute(
                protocols=protocols,
                clients_indices=[0] * RetriesTest._n_protocols,
                resources=None,
                priorities=None,
                retries=RetriesTest._n_retries_per_protocol,
            )
        # Only one counter file should have the maximum attempts, but because submitted tasks
        # run asynchronously in a potentially shuffled order, we must find which one
        counter_files = glob.glob(
            os.path.join(output_path, RetriesTest._counter_filename_template.format(protocol_number="*", task="*"))
        )
        task_total_attempts_dict = {}
        for counter_file in counter_files:
            basename_split = os.path.basename(counter_file).split("__")
            protocol_number = int(list(filter(lambda x: x.startswith("protocol_number"), basename_split))[0].split("-")[-1])
            self.assertEqual(protocol_number, 0, msg="Only protocol number 0 should have been run!")
            task = int(list(filter(lambda x: x.startswith("task"), basename_split))[0].split("-")[-1])
            self.assertIn(task, tuple(range(RetriesTest._n_tasks)), msg="Task is out of task range!")
            with open(counter_file, "r") as f:
                total_attempts = f.read().count(RetriesTest._count_marker)
            task_total_attempts_dict[task] = total_attempts
        task_with_max_attempts, max_attempts = sorted(task_total_attempts_dict.items(), key=lambda kv: kv[1], reverse=True)[0]
        self.assertEqual(max_attempts, RetriesTest._n_retries_per_protocol + 1, msg="Maximum number of task attempts failed.")
        for task, attempts in task_total_attempts_dict.items():
            if task != task_with_max_attempts:
                self.assertLess(
                    attempts,
                    max_attempts,
                    msg=(
                        f"Task {task} with {attempts} attempts must have less attempts "
                        f"than task {task_with_max_attempts} with maximum attempts {max_attempts}."
                    ),
                )
        print("Unit test with protocols raising intentional errors passed successfully!")

    def test_retries_succeed_on_last_retry(self):
        """Test retries of protocols that succeed on last retry."""
        protocols = [RetriesTest.protocol_succeeds_on_last_retry] * RetriesTest._n_protocols
        output_path = os.path.join(self.workdir.name, "outputs_protocols_succeed_on_last_retry")
        PyRosettaCluster(
            **self.default_instance_kwargs,
            **{"output_path": output_path},
        ).distribute(
            protocols=protocols,
            clients_indices=[0] * RetriesTest._n_protocols,
            resources=None,
            priorities=None,
            retries=[RetriesTest._n_retries_per_protocol] * RetriesTest._n_protocols,
        )
        for protocol_number in range(RetriesTest._n_protocols):
            for task in range(RetriesTest._n_tasks):
                counter_file = os.path.join(
                    output_path,
                    RetriesTest._counter_filename_template.format(protocol_number=protocol_number, task=task),
                )
                self.assertTrue(os.path.isfile(counter_file), msg=f"Retries counter file does not exist: {counter_file}")
                with open(counter_file, "r") as f:
                    total_attempts = f.read().count(RetriesTest._count_marker)
                self.assertEqual(total_attempts, RetriesTest._n_retries_per_protocol + 1, msg="Number of task attempts failed.")
        print("Unit test with protocols succeeding only on last retry passed successfully!")

    def test_no_retries(self):
        """Test no retries of a failed protocol."""
        protocols = [RetriesTest.protocol_with_intentional_error] * RetriesTest._n_protocols
        output_path = os.path.join(self.workdir.name, "outputs_protocols_no_retries")
        prc = PyRosettaCluster(
            **self.default_instance_kwargs,
            **{"output_path": output_path},
        )
        with self.assertRaises(WorkerError):
            prc.distribute(
                protocols=protocols,
                clients_indices=[0] * RetriesTest._n_protocols,
                resources=None,
                priorities=None,
                retries=None, # Automatically uses the `dask.distributed` default of 0
            )
        # Because submitted tasks run asynchronously in a potentially shuffled order,
        # only one task should have been run
        tasks_run_counter = 0
        for protocol_number in range(RetriesTest._n_protocols):
            for task in range(RetriesTest._n_tasks):
                counter_file = os.path.join(
                    output_path,
                    RetriesTest._counter_filename_template.format(protocol_number=protocol_number, task=task),
                )
                if os.path.isfile(counter_file):
                    with open(counter_file, "r") as f:
                        total_attempts = f.read().count(RetriesTest._count_marker)
                    self.assertEqual(total_attempts, 1, msg="Number of task attempts failed.")
                    tasks_run_counter += 1
        self.assertEqual(tasks_run_counter, 1, msg=f"Unexpected number of tasks run: {tasks_run_counter}")
        print("Unit test with no retries passed successfully!")

    def test_retries_api(self):
        """Test retries API."""
        protocols = [RetriesTest.my_protocol] * RetriesTest._n_protocols
        clients_indices = [0] * RetriesTest._n_protocols
        output_path = os.path.join(self.workdir.name, "outputs_protocols_retries_api")
        produce(
            **self.default_instance_kwargs,
            **dict(
                output_path=f"{output_path}_1",
                protocols=protocols,
                clients_indices=clients_indices,
                resources=None,
                priorities=None,
                retries=tuple(range(0, RetriesTest._n_protocols * 10, 10)), # Test different values in tuple
            )
        )
        run(
            **self.default_instance_kwargs,
            **dict(
                output_path=f"{output_path}_2",
                protocols=protocols,
                clients_indices=clients_indices,
                resources=None,
                priorities=None,
                retries=list(reversed(range(RetriesTest._n_protocols))), # Test different values in list
            )
        )
        prc_iterator = iterate(
            **self.default_instance_kwargs,
            **dict(
                output_path=f"{output_path}_3",
                protocols=protocols,
                clients_indices=clients_indices,
                resources=None,
                priorities=None,
                retries=set(range(RetriesTest._n_protocols)), # Test validation error with set
            )
        )
        with self.assertRaises(ValueError) as ex:
            _ = list(prc_iterator)
        self.assertIn(
            "The `retries` keyword argument parameter must be of type `int`, `list`, or `tuple`.",
            str(ex.exception),
            msg="Incorrect error message."
        )
        prc_generator = PyRosettaCluster(
            **self.default_instance_kwargs,
            output_path=f"{output_path}_4",
        ).generate(
            protocols=protocols,
            clients_indices=clients_indices,
            resources=None,
            priorities=None,
            retries=list(range(RetriesTest._n_protocols + 1)), # Test validation error with different size
        )
        with self.assertRaises(ValueError) as ex:
            _ = list(prc_generator)
        self.assertIn(
            "The `retries` keyword argument parameter must have the same length as the number of user-defined PyRosetta protocols!",
            str(ex.exception),
            msg="Incorrect error message."
        )
        prc = PyRosettaCluster(
            **self.default_instance_kwargs,
            output_path=f"{output_path}_5",
        )
        with self.assertRaises(ValueError) as ex:
            prc.distribute(
                protocols=protocols,
                clients_indices=clients_indices,
                resources=None,
                priorities=None,
                retries=-1, # Test validation error with negative integer
            )
        self.assertIn(
            "If the `retries` keyword argument parameter is of type `int`, it must be greater than or equal to 0.",
            str(ex.exception),
            msg="Incorrect error message."
        )
        with self.assertRaises(ValueError) as ex:
            produce(
                **self.default_instance_kwargs,
                **dict(
                    output_path=f"{output_path}_6",
                    protocols=protocols,
                    clients_indices=clients_indices,
                    resources=None,
                    priorities=None,
                    retries=1.0, # Test validation error with float
                )
            )
        self.assertIn(
            "The `retries` keyword argument parameter must be of type `int`, `list`, or `tuple`.",
            str(ex.exception),
            msg="Incorrect error message."
        )
        print("Unit test for retries API passed successfully!")


@unittest.skipIf(
    sys.version_info[:2] < (3, 8),
    "Unit test implements positional-only argument separators in function signatures, which is supported by Python-3.8+ only."
)
class TaskKeysTest(unittest.TestCase):
    """Smoke tests for user-defined task dictionary keys and user-defined PyRosetta protocol input signatures."""
    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        cls.input_packed_pose = io.pose_from_sequence("ARG/SEP")
        cls.workdir = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    @staticmethod
    def positional_arbitrary_parameter_protocol(baz, /, **kwargs):
        """Test Case #1: A user-defined PyRosetta protocol with a positional-only parameter named 'baz'."""
        return baz, kwargs

    @staticmethod
    def positional_or_keyword_arbitrary_parameter_protocol(baz, **kwargs):
        """Test Case #2: A user-defined PyRosetta protocol with a positional-or-keyword parameter named 'baz'."""
        return baz, kwargs

    @staticmethod
    def positional_custom_parameter_protocol(foo, /, **kwargs):
        """Test Case #3: A user-defined PyRosetta protocol with a positional-only parameter named 'foo'."""
        return foo, kwargs

    @staticmethod
    def positional_or_keyword_custom_parameter_protocol(foo, **kwargs):
        """Test Case #4: A user-defined PyRosetta protocol with a positional-or-keyword parameter named 'foo'."""
        return foo, kwargs

    @staticmethod
    def positional_parameter_protocol(packed_pose, /, **kwargs):
        """Test Case #5: A user-defined PyRosetta protocol with a positional-only parameter named 'packed_pose'."""
        return packed_pose, kwargs

    @staticmethod
    def positional_or_keyword_parameter_protocol(packed_pose, **kwargs):
        """Test Case #6: A user-defined PyRosetta protocol with a positional-or-keyword parameter named 'packed_pose'."""
        return packed_pose, kwargs

    def test_task_keys(self):
        """Test compatibility between user-defined task dictionary keys and user-defined PyRosetta protocol input signatures."""
        task_keys = set()
        for func in (run_protocol, user_protocol):
            sig = inspect.signature(func)
            for name, param in sig.parameters.items():
                if str(param.kind) != "VAR_KEYWORD": # Add all except variadic keyword parameters
                    task_keys.add(name)
        task_keys.update({"foo", "bar"}) # Add custom task keys
        instance_kwargs = dict(
            tasks=dict(map(reversed, enumerate(task_keys))),
            input_packed_pose=self.input_packed_pose,
            scheduler=None,
            cores=None,
            processes=None,
            memory=None,
            min_workers=1,
            max_workers=1,
            nstruct=1,
            project_name="PyRosettaCluster_Tests",
            simulation_name=uuid.uuid4().hex,
            scratch_dir=os.path.join(self.workdir.name, "scratch"),
            sha1=None,
        )
        output_path = os.path.join(self.workdir.name, "test_task_keys")
        _sep = "*" * 60

        # --- Test case #1 ---
        # Test a user-defined PyRosetta protocol defining a positional-only parameter named 'baz':
        # The user-defined task dictionary does not contain the 'baz' key, and the user-defined PyRosetta
        # protocol defines a positional-only parameter named 'baz', so upon calling the user-defined PyRosetta protocol,
        # PyRosettaCluster passes the `packed_pose` variable as a positional argument to the 'baz' parameter, and the
        # 'packed_pose' keyword argument from the unpacked user-defined task dictionary is collected by the variadic
        # keyword parameter, which avoids a `TypeError` since it's valid Python syntax.
        protocol = TaskKeysTest.positional_arbitrary_parameter_protocol
        parameter_name = "baz"
        self.assertNotIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_1")

        # --- Test case #2 ---
        # Test a user-defined PyRosetta protocol defining a positional-or-keyword parameter named 'baz':
        # The user-defined task dictionary does not contain the 'baz' key, and the user-defined PyRosetta
        # protocol defines a positional-or-keyword parameter named 'baz', so upon calling the user-defined PyRosetta protocol,
        # PyRosettaCluster passes the `packed_pose` variable as a positional argument to the 'baz' parameter, and the
        # 'packed_pose' keyword argument from the unpacked user-defined task dictionary is collected by the variadic
        # keyword parameter, which avoids a `TypeError` since it's valid Python syntax.
        protocol = TaskKeysTest.positional_or_keyword_arbitrary_parameter_protocol
        parameter_name = "baz"
        self.assertNotIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_2")

        # --- Test case #3 ---
        # Test a user-defined PyRosetta protocol defining a positional-only parameter named 'foo':
        # Although the user-defined task dictionary contains the 'foo' key, because the user-defined PyRosetta
        # protocol defines a positional-only parameter named 'foo', upon calling the user-defined PyRosetta protocol,
        # PyRosettaCluster passes the `packed_pose` variable as a positional argument to the 'foo' parameter, and the
        # 'foo' keyword argument from the unpacked user-defined task dictionary is collected by the variadic
        # keyword parameter, which avoids a `TypeError` since it's valid Python syntax.
        protocol = TaskKeysTest.positional_custom_parameter_protocol
        parameter_name = "foo"
        self.assertIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_3")

        # --- Test case #4 ---
        # Test a user-defined PyRosetta protocol defining a positional-or-keyword parameter named 'foo':
        protocol = TaskKeysTest.positional_or_keyword_custom_parameter_protocol
        parameter_name = "foo"
        self.assertIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        print(f"{_sep} Begin testing expected `TypeError` in billiard subprocess {_sep}")
        with self.assertRaises(WorkerError):
            # Catch exception in billiard subprocess raised due to the user-defined task dictionary containing the
            # 'foo' key and the user-defined PyRosetta protocol defining the positional-or-keyword parameter 'foo'.
            # When calling the user-defined PyRosetta protocol, PyRosettaCluster passes the `packed_pose` variable as
            # a positional argument to the 'foo' parameter, and the user-defined task dictionary is unpacked with the
            # 'foo' keyword argument name, causing the exception:
            #     TypeError: TaskKeysTest.positional_or_keyword_custom_parameter_protocol() got multiple values for argument 'foo'
            # However, this is a user error. To fix it, the user must either: (i) remove the 'foo' key from the user-defined
            # task dictionary; (ii) add a positional-only argument separator to the user-defined PyRosetta protocol input
            # signature to make the 'foo' parameter positional-only (as in Test Case #3); or (iii) change the positional-or-keyword
            # parameter name in the user-defined PyRosetta protocol input signature so that it doesn't clash with the 'foo' key
            # from the user-defined task dictionary (as in Test Case #2).
            run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_4")
        print(f"{_sep} End testing expected `TypeError` in billiard subprocess {_sep}")

        # --- Test case #5 ---
        # Test a user-defined PyRosetta protocol defining a positional-only parameter named 'packed_pose':
        # Although the user-defined task dictionary contains the 'packed_pose' key, because the user-defined PyRosetta
        # protocol defines a positional-only parameter named 'packed_pose', upon calling the user-defined PyRosetta protocol,
        # PyRosettaCluster passes the `packed_pose` variable as a positional argument to the 'packed_pose' parameter, and the
        # 'packed_pose' keyword argument from the unpacked user-defined task dictionary is collected by the variadic
        # keyword parameter, which avoids a `TypeError` since it's valid Python syntax.
        protocol = TaskKeysTest.positional_parameter_protocol
        parameter_name = "packed_pose"
        self.assertIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_5")

        # --- Test case #6 ---
        # Test a user-defined PyRosetta protocol defining a positional-or-keyword parameter named 'packed_pose':
        protocol = TaskKeysTest.positional_or_keyword_parameter_protocol
        parameter_name = "packed_pose"
        self.assertIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        print(f"{_sep} Begin testing expected `TypeError` in billiard subprocess {_sep}")
        with self.assertRaises(WorkerError):
            # Catch exception in billiard subprocess raised due to the user-defined task dictionary containing the
            # 'packed_pose' key and the user-defined PyRosetta protocol defining the positional-or-keyword parameter 'packed_pose'.
            # When calling the user-defined PyRosetta protocol, PyRosettaCluster passes the `packed_pose` variable as a positional
            # argument to the 'packed_pose' parameter, and the user-defined task dictionary is unpacked with the 'packed_pose'
            # keyword argument name, causing the exception:
            #     TypeError: TaskKeysTest.positional_or_keyword_parameter_protocol() got multiple values for argument 'packed_pose'
            # However, this is a user error. To fix it, the user must either: (i) remove the 'packed_pose' key from the
            # user-defined task dictionary; (ii) add a positional-only argument separator to the user-defined PyRosetta protocol
            # input signature to make the 'packed_pose' parameter positional-only (as in Test Case #5); or (iii) change the
            # positional-or-keyword parameter name in the user-defined PyRosetta protocol input signature so that it doesn't clash
            # with the 'packed_pose' key from the user-defined task dictionary (as in Test Case #2).
            run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_6")
        print(f"{_sep} End testing expected `TypeError` in billiard subprocess {_sep}")


class WorkerPreemptionTest(unittest.TestCase):
    """Smoke tests for preempted compute resources with PyRosettaCluster."""
    # `WorkerPreemptionTest._allowed_failures` must be a very large number for `WorkerPreemptionTest` unit tests,
    # otherwise task chains are terminated on client upon catching `KilledWorker` exceptions and test logic will fail.
    _allowed_failures = 9999
    _n_workers = 4 # Must be >1 for `WorkerPreemptionTest`
    # `WorkerPreemptionTest._sleep_time` must be sufficient to allow Dask to have enough time to reschedule lost tasks on alive workers
    # before the next worker preemption. Empirically, the required sleep time depends on: the number of tasks in flight, the scheduler
    # heartbeat, PyRosettaCluster overhead time (e.g., serializing data), etc.
    _sleep_time = 3
    _n_tasks = 4
    _n_protocols = 2
    _sep = "*" * 60

    @classmethod
    def setUpClass(cls):
        random.seed(12345)
        # Increase resilience to worker preemption
        dask.config.set(
            {
                "distributed.scheduler.allowed-failures": WorkerPreemptionTest._allowed_failures,
                "distributed.scheduler.active-memory-manager.start": False, # Must be disabled for `WorkerPreemptionTest` test case
            }
        )
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        cls.workdir = tempfile.TemporaryDirectory()
        with warnings.catch_warnings():
            # Catch 'ResourceWarning: unclosed <socket.socket ...' from distributed/node.py:235
            # Catch 'UserWarning: Port 8787 is already in use' from distributed/node.py:240
            # Catch 'DeprecationWarning: `np.bool8` is a deprecated alias for `np.bool_`.  (Deprecated NumPy 1.24)' from bokeh/core/property/primitive.py:37
            # Catch 'DeprecationWarning: pkg_resources is deprecated as an API.' from jupyter_server_proxy/config.py:10
            warnings.simplefilter("ignore", category=ResourceWarning)
            warnings.simplefilter("ignore", category=UserWarning)
            warnings.simplefilter("ignore", category=DeprecationWarning)
            cls.cluster_1 = LocalCluster(
                n_workers=WorkerPreemptionTest._n_workers,
                threads_per_worker=1,
                memory_limit="8GB", # Must be sufficient enough to support potentially all tasks replicated on all workers (with ~1MB/task)
                dashboard_address=None,
                local_directory=cls.workdir.name,
            )
        cls.client_1 = Client(cls.cluster_1)
        cls.scorefile_name = "scores.json"
        cls.dry_run = False # Must be `False` for `WorkerPreemptionTest` test case to read the output scorefile
        cls.save_all = True # Must be `True` for `WorkerPreemptionTest` test case to preempt Dask worker processes mid-trajectory
        cls.simulation_records_in_scorefile = True # Must be `True` for `WorkerPreemptionTest` test case to read the output scorefile
        cls.default_instance_kwargs = dict(
            tasks=WorkerPreemptionTest.create_tasks,
            input_packed_pose=io.pose_from_sequence("ACDEFGHIKLMNPQRSTVWY" * 50),
            seeds=None,
            decoy_ids=None,
            client=None,
            clients=[cls.client_1],
            scheduler=None,
            scratch_dir=cls.workdir.name,
            cores=None,
            processes=None,
            memory=None,
            min_workers=1,
            max_workers=1,
            nstruct=1,
            dashboard_address=None,
            compressed=True,
            logging_level="INFO",
            scorefile_name=cls.scorefile_name,
            project_name="PyRosettaCluster_Tests",
            simulation_name=uuid.uuid4().hex,
            environment=None,
            simulation_records_in_scorefile=cls.simulation_records_in_scorefile,
            decoy_dir_name="decoys",
            logs_dir_name="logs",
            ignore_errors=False,
            timeout=0.5,
            max_delay_time=0.5,
            sha1=None,
            dry_run=cls.dry_run,
            save_all=cls.save_all,
            system_info=None,
            pyrosetta_build=None,
            filter_results=True,
            norm_task_options=None,
            output_init_file=None,
            output_decoy_types=[".pdb"],
            output_scorefile_types=[".json"],
            security=False,
            max_nonce=99999,
        )
        print(
            WorkerPreemptionTest._sep,
            "Start testing worker preemption with expected 'Restarting worker'/'Removing worker' warnings",
            WorkerPreemptionTest._sep,
        )

    @classmethod
    def tearDownClass(cls):
        cls.client_1.close()
        cls.cluster_1.close()
        cls.workdir.cleanup()
        print(
            WorkerPreemptionTest._sep,
            "End testing worker preemption with expected 'Restarting worker'/'Removing worker' warnings",
            WorkerPreemptionTest._sep,
        )

    @staticmethod
    def create_tasks():
        for i in range(WorkerPreemptionTest._n_tasks):
            yield {
                "extra_options": "-ex1 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
                "task": i,
            }

    @staticmethod
    def my_protocol(packed_pose, **kwargs):
        return packed_pose, kwargs

    def get_current_worker_pids_map(self):
        return {k: v.pid for k, v in self.client_1.cluster.scheduler.workers.items()}

    def preempt_a_worker(self, non_preemptible_worker_pids_map, max_attempts=20, verbose=True):
        """
        Randomly select and terminate a worker process and its billiard subprocesses
        to simulate compute resource preemption. Method adapted from:
        https://examples.dask.org/resilience.html#Suddenly-shutting-down-workers
        """
        def _terminate_process_tree(pid, sig=signal.SIGTERM):
            try:
                parent = psutil.Process(pid)
            except psutil.NoSuchProcess:
                return False
            for child in parent.children(recursive=True):
                try:
                    child.send_signal(sig)
                except psutil.NoSuchProcess:
                    pass
            try:
                parent.send_signal(sig)
                return True
            except psutil.NoSuchProcess:
                return False

        _attempts = 0
        while _attempts < max_attempts:
            _attempts += 1
            current_worker_pids_map = self.get_current_worker_pids_map()
            if verbose:
                print("Current workers and process IDs:")
                pprint(current_worker_pids_map)
            preemptible_worker_pids_map = toolz.dicttoolz.keyfilter(
                lambda k: k not in non_preemptible_worker_pids_map,
                current_worker_pids_map
            )
            if not preemptible_worker_pids_map:
                if verbose:
                    print(f"Warning: attempt {_attempts} failed because no preemptible Dask workers are available. Retrying...")
                time.sleep(0.5)
                continue
            _worker = random.choice(list(preemptible_worker_pids_map.keys()))
            _worker_pid = preemptible_worker_pids_map[_worker]
            _terminated = _terminate_process_tree(_worker_pid)
            if _terminated:
                if verbose:
                    print(f"Terminated worker process ID: {_worker_pid}")
                    print(f"Terminated worker: {_worker}")
                time.sleep(WorkerPreemptionTest._sleep_time)
                break
            else:
                if verbose:
                    print(f"Warning: attempt {_attempts} failed to terminate Dask worker '{_worker}' with process ID {_worker_pid}. Retrying...")
            time.sleep(0.5)
        else:
            raise RuntimeError(f"No Dask workers available to preempt after {_attempts} attempts: {current_worker_pids_map}")

    def simulate_worker_preemption(
        self,
        max_task_replicas=0,
        task_registry=None,
        n_non_preemptible_workers=2, # Must be < `WorkerPreemptionTest._n_workers`
        verbose=False,
    ):
        self.assertEqual(dask.config.get("distributed.scheduler.allowed-failures"), WorkerPreemptionTest._allowed_failures)
        self.assertFalse(dask.config.get("distributed.scheduler.active-memory-manager.start"))
        # Mark some worker processes as non-preemptible to allow task chains to eventually run to completion
        self.assertGreater(WorkerPreemptionTest._n_workers, 1)
        if not task_registry:
            self.assertGreater(n_non_preemptible_workers, 1)
            self.assertLess(n_non_preemptible_workers, WorkerPreemptionTest._n_workers)
        current_worker_pids_map = self.get_current_worker_pids_map()
        non_preemptible_worker_pids_map = dict(
            toolz.itertoolz.take(
                n_non_preemptible_workers,
                current_worker_pids_map.items(),
            )
        )
        if verbose:
            print(f"Non-preemptible workers and process IDs: {non_preemptible_worker_pids_map}")
        # Setup simulation
        output_path = os.path.join(
            self.workdir.name,
            f"outputs_worker_preemption_{max_task_replicas}_{task_registry}_{n_non_preemptible_workers}",
        )
        prc_iterator = iterate(
            **self.default_instance_kwargs,
            output_path=output_path,
            max_task_replicas=max_task_replicas,
            task_registry=task_registry,
        )
        prc = PyRosettaCluster(
            **self.default_instance_kwargs,
            output_path=output_path,
            max_task_replicas=max_task_replicas,
            task_registry=task_registry,
        )
        self.assertFalse(getattr(prc, "dry_run", None))
        self.assertTrue(getattr(prc, "save_all", None))
        self.assertTrue(getattr(prc, "simulation_records_in_scorefile", None))
        if task_registry is None:
            self.assertEqual(prc.registry, None)
        elif task_registry == "disk":
            self.assertIsInstance(prc.registry, DiskTaskRegistry)
        elif task_registry == "memory":
            self.assertIsInstance(prc.registry, MemoryTaskRegistry)
        prc_iterator = prc.generate(
            protocols=[WorkerPreemptionTest.my_protocol] * WorkerPreemptionTest._n_protocols,
            clients_indices=[0] * WorkerPreemptionTest._n_protocols,
            resources=None,
            priorities=None,
            retries=None,
        )
        # Start the simulation and keep shutting down workers while it is running using the
        # `PyRosettaCluster.generate` method to intervene periodically throughout the simulation
        self.assertGreaterEqual(WorkerPreemptionTest._sleep_time, 3)
        for _packed_pose, _kwargs in prc_iterator:
            if verbose:
                print("Yielded decoy IDs:", _kwargs.get("PyRosettaCluster_decoy_ids"))
                print("Who has:")
                pprint(self.client_1.who_has())
                print("Processing:")
                pprint(self.client_1.processing())
                if task_registry:
                    print(f"Current size of task registry (MB): {prc.registry.total_size() / 1e6:.6f}")
                    print(f"Current number of task registry entries: {len(prc.registry)}")
                    print("Current task registry keys:")
                    pprint(list(prc.registry))
            self.preempt_a_worker(non_preemptible_worker_pids_map, verbose=verbose)
        # Assert that all task chains completed
        scorefile = os.path.join(output_path, self.scorefile_name)
        with open(scorefile, "r") as f:
            _n_results = sum(1 for line in f)
        _n_expected_results = WorkerPreemptionTest._n_tasks * WorkerPreemptionTest._n_protocols
        self.assertEqual(
            _n_results,
            _n_expected_results,
            msg=f"Task chains did not run to completion with {n_non_preemptible_workers} non-preemptible workers.",
        )
        # Assert that full simulation records do not save the 'max_task_replicas' or 'task_registry' instance attributes
        with open(scorefile, "r") as f:
            for line in f:
                record = json.loads(line)
                for entry in ("instance", "metadata", "scores"):
                    self.assertNotIn("max_task_replicas", record[entry])
                    self.assertNotIn("task_registry", record[entry])
                break

    def test_disk_task_registry(self):
        self.simulate_worker_preemption(
            max_task_replicas=0,
            task_registry="disk",
            n_non_preemptible_workers=0,
        )

    def test_memory_task_registry(self):
        self.simulate_worker_preemption(
            max_task_replicas=0,
            task_registry="memory",
            n_non_preemptible_workers=0,
        )

    def test_disk_max_task_replicas_all(self):
        self.simulate_worker_preemption(
            max_task_replicas=None,
            task_registry="disk", # Fallback, otherwise unit test might rarely fail stochastically
            n_non_preemptible_workers=WorkerPreemptionTest._n_workers // 2,
        )

    def test_disk_max_task_replicas_int(self):
        self.simulate_worker_preemption(
            max_task_replicas=WorkerPreemptionTest._n_workers // 2,
            task_registry="memory", # Fallback, otherwise unit test might fail stochastically
            n_non_preemptible_workers=WorkerPreemptionTest._n_workers // 2,
        )
