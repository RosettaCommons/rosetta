# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"
__email__ = "klima.jason@gmail.com"

try:
    import psutil
    from dask.distributed import Adaptive, Client, LocalCluster
    from dask_jobqueue import SGECluster, SLURMCluster
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

from typing import (
    Generic,
    Optional,
    Tuple,
    TypeVar,
)


AdaptiveType = TypeVar("AdaptiveType", bound=Adaptive)
ClientType = TypeVar("ClientType", bound=Client)
ClusterType = TypeVar("ClusterType", LocalCluster, SGECluster, SLURMCluster)
G = TypeVar("G")


class SchedulerManager(Generic[G]):
    def _get_cluster(self) -> ClusterType:
        """Given user input arguments, return the requested cluster instance."""

        if not self.scheduler:
            _cpu_count = psutil.cpu_count()
            _n_workers = (
                _cpu_count if (_cpu_count < self.max_workers) else self.max_workers
            )
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

    def _setup_client_cluster_adaptive(
        self,
    ) -> Tuple[
        ClientType, Optional[ClusterType], Optional[AdaptiveType],
    ]:
        """
        Given user input arguments, return the requested client, cluster,
        and adaptive instance.
        """
        if self.client:
            client = self.client
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

        return client, cluster, adaptive

    def _maybe_adapt(self, adaptive: Optional[AdaptiveType]) -> None:
        """Adjust max_workers."""

        if (
            not self.client
            and self.scheduler
            and (self.max_workers >= 1000)
            and adaptive
        ):
            adaptive.maximum = self.tasks_size

    def _maybe_teardown(
        self, client: ClientType, cluster: Optional[ClusterType],
    ) -> None:
        """Teardown client and cluster."""

        logging.info("PyRosettaCluster simulation complete!")
        if not self.client and cluster:
            cluster.scale(0)
            client.close()
            cluster.close()
