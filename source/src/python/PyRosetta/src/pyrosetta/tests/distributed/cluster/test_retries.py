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
import time
import unittest
import uuid
import warnings

try:
    from dask.distributed import Client, LocalCluster
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_retries' requires the "
        + "third-party package 'dask' as a dependency!\n"
        + "Please install this package into your virtual environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n",
        flush=True,
    )
    raise

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    iterate,
    produce,
    run,
)
from pyrosetta.distributed.cluster.exceptions import WorkerError
from pyrosetta.tests.distributed.cluster.unittest_utils import IntentionalError


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
        print(f"{RetriesTest._sep} Begin testing PyRosettaCluster().distribute(retries=...) {RetriesTest._sep}", flush=True)

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
        print(f"{RetriesTest._sep} End testing PyRosettaCluster().distribute(retries=...) {RetriesTest._sep}", flush=True)

    def tearDown(self):
        sys.stdout.flush()

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
            print(f"Protocol `protocol_succeeds_on_last_retry` succeeding on last attempt ({total_attempts}/{max_attempts})!", flush=True)
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
        print("Unit test with protocols raising intentional errors passed successfully!", flush=True)

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
        print("Unit test with protocols succeeding only on last retry passed successfully!", flush=True)

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
        print("Unit test with no retries passed successfully!", flush=True)

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
            "The `retries` keyword argument value must be of type `int`, `list`, or `tuple`.",
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
            "The `retries` keyword argument value must have the same length as the number of user-defined PyRosetta protocols!",
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
            "If the `retries` keyword argument value is of type `int`, it must be greater than or equal to 0.",
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
            "The `retries` keyword argument value must be of type `int`, `list`, or `tuple`.",
            str(ex.exception),
            msg="Incorrect error message."
        )
        print("Unit test for retries API passed successfully!", flush=True)
