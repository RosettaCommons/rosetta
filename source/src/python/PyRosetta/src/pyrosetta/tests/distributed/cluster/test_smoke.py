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
import tempfile
import unittest
import uuid

try:
    import cryptography
    has_cryptography = True
except ImportError as ex:
    has_cryptography = False

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    produce,
    run,
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
            security = pyrosetta.distributed.cluster.generate_dask_tls_security(
                os.path.join(workdir, "security_test")
            )
            print(
                "Successfully ran `pyrosetta.distributed.cluster.generate_dask_tls_security()` "
                + f"to generate a dask `Security` object: {security}",
                flush=True,
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
                print("Using the installed 'cryptography' package to generate a dask `Security.temporary()` object...", flush=True)
                run(**instance_kwargs)
            else:
                with self.assertRaises(ImportError) as ex:
                    run(**instance_kwargs)
                print(f"Successfully caught exception: {ex.exception}", flush=True)

    def test_invalid_tasks(self):
        """Smoke test for catching invalid input tasks in PyRosettaCluster."""
        pyrosetta.distributed.init("-run:constant_seed 1 -multithreading:total_threads 1")
        with tempfile.TemporaryDirectory() as workdir:
            _invalid_run_options_values = (
                ("constant_seed", "1"),
                ("jran", "12345"),
                ("use_time_as_seed", "1"),
                ("rng_seed_device", "/dev/urandom"),
                ("seed_offset", "2"),
                ("rng", "mt19937"),
            )
            _run_prefixes = ("-run::", "-run:", "-")
            _options_keys = ("options", "extra_options")
            _invalid_options_values = [
                f"{prefix}{k} {v}"
                for prefix in _run_prefixes
                for k, v in _invalid_run_options_values
            ]
            _invalid_options_values += [
                {k: v} for k, v in map(lambda val: val.split(), _invalid_options_values)
            ]
            _invalid_options_values += [True, 123, 123.456, complex(1, 2), b"foo"]
            _invalid_tasks = [
                {_option: _value}
                for _option in _options_keys
                for _value in _invalid_options_values
            ]
            _invalid_tasks += [
                {"PyRosettaCluster_foo": "bar"},
                {True: "bar"},
                {0: "bar"},
                {b"Bytes": "bar"},
            ]
            for _invalid_task in _invalid_tasks:
                for _option in _options_keys:
                    if _option not in _invalid_task:
                        _invalid_task[_option] = ""
                with self.assertRaises(ValueError) as ex:
                    PyRosettaCluster(
                        tasks=_invalid_task,
                        output_path=os.path.join(workdir, f"outputs_invalid_task_{uuid.uuid4().hex}"),
                        scratch_dir=workdir,
                        sha1=None,
                    )
                print(f"Successfully caught exception: {ex.exception}", flush=True)

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
        print(f"{_sep} Begin testing PyRosettaCluster(ignore_errors=...) {_sep}", flush=True)

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

        print(f"{_sep} End testing PyRosettaCluster(ignore_errors=...) {_sep}", flush=True)
