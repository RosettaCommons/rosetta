# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import psutil
    from dask.distributed import Adaptive, Client, LocalCluster
    from dask_jobqueue import SGECluster, SLURMCluster
    from distributed.protocol import dask_serialize, dask_deserialize
    from distributed.protocol.serialize import register_serialization_family
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.utilities' requires the "
        + "third-party packages 'dask', 'dask-jobqueue', and 'psutil' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
        + "https://pypi.org/project/dask-jobqueue/\n"
        + "https://pypi.org/project/psutil/\n"
    )
    raise

import logging
import os
import warnings

from typing import (
    Dict,
    Generic,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
)


AdaptiveType = TypeVar("AdaptiveType", bound=Adaptive)
ClientType = TypeVar("ClientType", bound=Client)
ClusterType = TypeVar("ClusterType", LocalCluster, SGECluster, SLURMCluster)
G = TypeVar("G")


class SchedulerManager(Generic[G]):
    def _setup_clients_dict(self) -> Union[Dict[int, ClientType], NoReturn]:
        if all(x is None for x in (self.client, self.clients)):
            return {}
        elif isinstance(self.client, Client) and self.clients is None:
            return {0: self.client}
        elif isinstance(self.clients, (tuple, list)) and self.client is None:
            return dict(enumerate(self.clients, start=0))
        else:
            raise ValueError(
                "The PyRosettaCluster `client` and `clients` attribute parameters may not both be set. Received:\n" 
                + f"PyRosettaCluster().client: {self.client}\n"
                + f"PyRosettaCluster().clients: {self.clients}\n"
            )

    def _get_cluster(self) -> ClusterType:
        """Given user input arguments, return the requested cluster instance."""

        if not self.scheduler:
            _cpu_count = psutil.cpu_count()
            _n_workers = (
                _cpu_count if (_cpu_count < self.max_workers) else self.max_workers
            )
            with warnings.catch_warnings():
                # Catch 'ResourceWarning: unclosed <socket.socket ...' from distributed/node.py:235
                # Catch 'UserWarning: Port 8787 is already in use' from distributed/node.py:240
                # Catch 'DeprecationWarning: `np.bool8` is a deprecated alias for `np.bool_`.  (Deprecated NumPy 1.24)' from bokeh/core/property/primitive.py:37
                warnings.simplefilter("ignore", category=ResourceWarning)
                warnings.simplefilter("default", category=UserWarning)
                warnings.simplefilter("ignore", category=DeprecationWarning)
                cluster = LocalCluster(
                    n_workers=_n_workers,
                    threads_per_worker=1,
                    dashboard_address=self.dashboard_address,
                    local_directory=self.scratch_dir,
                )
        else:
            if self.scheduler == "sge":
                cluster_func = SGECluster
                log_files = self.logs_path
            elif self.scheduler == "slurm":
                cluster_func = SLURMCluster
                log_files = os.path.join(self.logs_path, "dask-worker.o%j")
            cluster = cluster_func(
                cores=self.cores,
                processes=self.processes,
                memory=self.memory,
                local_directory=self.scratch_dir,
                job_extra=[f"-o {log_files}",],
                walltime="99999:0:0",
                death_timeout=9999,
                dashboard_address=self.dashboard_address,
            )
        logging.info(f"Dashboard link: {cluster.dashboard_link}")

        return cluster

    def _setup_clients_cluster_adaptive(
        self,
    ) -> Tuple[
        Dict[int, ClientType], Optional[ClusterType], Optional[AdaptiveType],
    ]:
        """
        Given user input arguments, return the requested client, cluster,
        and adaptive instance.
        """
        if self.clients_dict:
            clients = self.clients_dict
            cluster = None
            adaptive = None
        else:
            cluster = self._get_cluster()
            client = Client(cluster)
            if self.scheduler:
                adaptive = cluster.adapt(
                    minimum=self.min_workers, maximum=self.max_workers, interval="5s",
                )
                cluster._adaptive.adapt()
            else:
                adaptive = None
            clients = {0: client}

        return clients, cluster, adaptive

    def _maybe_adapt(self, adaptive: Optional[AdaptiveType]) -> None:
        """Adjust max_workers."""

        if (
            not self.clients_dict
            and self.scheduler
            and (self.max_workers >= 1000)
            and adaptive
        ):
            adaptive.maximum = self.tasks_size

    def _maybe_teardown(
        self, clients: Dict[int, ClientType], cluster: Optional[ClusterType],
    ) -> None:
        """Teardown client and cluster."""

        logging.info("PyRosettaCluster simulation complete!")
        if not self.clients_dict and cluster:
            cluster.scale(0)
            clients[0].close()
            cluster.close()
