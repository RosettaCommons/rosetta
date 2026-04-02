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

import json
import os
import psutil
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import random
import signal
import sys
import tempfile
import time
import unittest
import uuid
import warnings

try:
    import dask
    import toolz
    from dask.distributed import Client, LocalCluster
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_worker_preemption' requires the "
        + "third-party packages 'dask' and 'toolz' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
        + "https://pypi.org/project/toolz/\n",
        flush=True,
    )
    raise

from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    iterate,
)
from pyrosetta.distributed.cluster.task_registry import DiskTaskRegistry, MemoryTaskRegistry
from pyrosetta.utility import pprint_flush


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
            flush=True,
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
            flush=True,
        )

    def tearDown(self):
        sys.stdout.flush()

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

    def preempt_a_worker(self, non_preemptible_worker_pids_map, max_attempts=50, verbose=False):
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

        self.assertGreaterEqual(WorkerPreemptionTest._sleep_time, 3)
        preempted = False
        _attempts = 0
        while _attempts < max_attempts:
            _attempts += 1
            current_worker_pids_map = self.get_current_worker_pids_map()
            if verbose:
                print("Current workers and process IDs:", flush=True)
                pprint_flush(current_worker_pids_map)
            preemptible_worker_pids_map = toolz.dicttoolz.keyfilter(
                lambda k: k not in non_preemptible_worker_pids_map,
                current_worker_pids_map
            )
            if not preemptible_worker_pids_map:
                if verbose:
                    print(f"Preemption attempt {_attempts} has no preemptible Dask workers. Skipping preemption.", flush=True)
                preempted = False
                break
            _worker = random.choice(list(preemptible_worker_pids_map.keys()))
            _worker_pid = preemptible_worker_pids_map[_worker]
            _terminated = _terminate_process_tree(_worker_pid)
            if _terminated:
                if verbose:
                    print(f"Terminated worker process ID: {_worker_pid}", flush=True)
                    print(f"Terminated worker: {_worker}", flush=True)
                preempted = True
                break
            else:
                if verbose:
                    print(f"Warning: preemption attempt {_attempts} failed to terminate Dask worker '{_worker}' with process ID {_worker_pid}. Retrying...", flush=True)
                time.sleep(random.uniform(0.5, 3.0))
        else:
            raise RuntimeError(
                f"No Dask workers available to preempt after {_attempts} attempts!\n"
                + f"Current Dask workers and PIDs: {current_worker_pids_map}\n"
                + f"Preemptible Dask workers and PIDs: {preemptible_worker_pids_map}\n"
                + f"Non-Preemptible Dask workers and PIDs: {non_preemptible_worker_pids_map}\n"
            )
        if verbose:
            print(f"Sleeping for {WorkerPreemptionTest._sleep_time} seconds...", flush=True)
        time.sleep(WorkerPreemptionTest._sleep_time)

        return preempted

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
            print(f"Non-preemptible workers and process IDs: {non_preemptible_worker_pids_map}", flush=True)
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
        total_preempted_workers = 0
        for _packed_pose, _kwargs in prc_iterator:
            if verbose:
                print("Yielded decoy IDs:", _kwargs.get("PyRosettaCluster_decoy_ids"), flush=True)
                print("Who has:", flush=True)
                pprint_flush(self.client_1.who_has())
                print("Processing:", flush=True)
                pprint_flush(self.client_1.processing())
                if task_registry:
                    print(f"Current size of task registry (MB): {prc.registry.total_size() / 1e6:.6f}", flush=True)
                    print(f"Current number of task registry entries: {len(prc.registry)}", flush=True)
                    print("Current task registry keys:", flush=True)
                    pprint_flush(list(prc.registry))
            preempted = self.preempt_a_worker(non_preemptible_worker_pids_map, verbose=verbose)
            total_preempted_workers += int(preempted)
        self.assertGreater(
            total_preempted_workers,
            0,
            msg=f"Finished PyRosettaCluster simulation, but preempted {total_preempted_workers} Dask workers in total."
        )
        if verbose:
            print(f"Finished PyRosettaCluster simulation. Preempted {total_preempted_workers} Dask workers in total.", flush=True)
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
