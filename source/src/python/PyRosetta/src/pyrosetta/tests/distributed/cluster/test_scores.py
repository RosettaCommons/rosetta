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
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import sys
import tempfile
import unittest
import uuid

try:
    import pandas
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_scores' requires the "
        + "third-party package 'pandas' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/pandas/\n",
        flush=True,
    )
    raise

from pyrosetta.distributed.cluster import (
    get_scores_dict,
    reserve_scores,
    run,
)
from pyrosetta.distributed.cluster.exceptions import WorkerError
from pyrosetta.distributed.cluster.io import _is_pandas_object_pyarrow_backed


class ScoresTest(unittest.TestCase):
    """Test case for scoring `PackedPose` objects in PyRosettaCluster."""
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

    def tearDown(self):
        sys.stdout.flush()

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
        secure_packages = list(pyrosetta.secure_unpickle.get_secure_packages())
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
        print(f"{_sep} Begin testing expected `UnpickleSecurityError` in billiard subprocess {_sep}", flush=True)
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
        print(f"{_sep} End testing expected `UnpickleSecurityError` in billiard subprocess {_sep}", flush=True)
