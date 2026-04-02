"""
PyRosettaCluster test utilities using the `unittest` framework.
"""
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import logging
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import sys
import tempfile
import uuid
import warnings

try:
    from dask.distributed import Client, LocalCluster
except ImportError:
    print(
        "Importing 'pyrosetta.tests.distributed.cluster.unittest_utils' requires the "
        + "third-party package 'dask' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
    )
    raise

from pprint import pprint
from pyrosetta.exceptions import PyRosettaIsNotInitializedError


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


class IntentionalError(RuntimeError):
    def __init__(self, *args):
        super().__init__(*args)


def score_function_is_available(name):
    if not pyrosetta.rosetta.basic.was_init_called():
        raise PyRosettaIsNotInitializedError("PyRosetta must be initialized to locate the database.")
    if not isinstance(name, str):
        raise ValueError("Score function name must be a `str` object.")
    if not name.endswith(".wts"):
        name += ".wts"
    db = pyrosetta.rosetta.basic.database.full_name("scoring/weights")

    return os.path.isfile(os.path.join(db, name))
