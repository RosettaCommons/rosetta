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

import collections
import json
import numpy
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
        "Importing 'pyrosetta.tests.distributed.cluster.test_priorities' requires the "
        + "third-party package 'dask' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n",
        flush=True,
    )
    raise

from pyrosetta.distributed.cluster import produce


class PrioritiesTest(unittest.TestCase):
    """Test case for task priorities in PyRosettaCluster."""

    def tearDown(self):
        sys.stdout.flush()

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
                print(f"Successfully tested `priorities` keyword argument: {_msg}", flush=True)

            client_1.close()
            cluster_1.close()
