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
import sys
import tempfile
import unittest
import uuid
import warnings

try:
    from dask.distributed import Client, LocalCluster
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.test_clients' requires the "
        + "third-party package 'dask' as a dependency!\n"
        + "Please install this package into your virtual environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n",
        flush=True,
    )
    raise

from pyrosetta.distributed.cluster import produce


class MultipleClientsTest(unittest.TestCase):
    """Test case for multiple clients with PyRosettaCluster."""

    def tearDown(self):
        sys.stdout.flush()

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
