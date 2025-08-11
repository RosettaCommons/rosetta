# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    from distributed import Worker, WorkerPlugin
except ImportError:
    try:
        from distributed import Worker
        from distributed.diagnostics.plugin import WorkerPlugin
    except ImportError:
        print(
            "Importing 'pyrosetta.distributed.cluster.worker_plugins' requires the "
            + "third-party package 'distributed' as a dependency!\n"
            + "Please install the package into your python environment. "
            + "For installation instructions, visit:\n"
            + "https://pypi.org/project/distributed/\n"
        )
        raise

import logging

from contextlib import suppress
from typing import Optional

from pyrosetta.distributed.cluster.logging_handlers import MultiSocketHandler
from pyrosetta.distributed.cluster.logging_filters import (
    DefaultProtocolNameFilter,
    DefaultSocketAddressFilter,
    DefaultTaskIdFilter,
)


SOCKET_LOGGER_PLUGIN_NAME: str = "PyRosettaCluster_socket_logger_plugin"
WORKER_LOGGER_NAME: str = "PyRosettaCluster_dask_worker"


class SocketLoggerPlugin(WorkerPlugin):
    """Install a `MultiSocketHandler` logging handler on a dask worker logger."""
    def __init__(self, logging_level: str, maxsize: int = 128):
        self.logging_level: str = logging_level
        self.maxsize: int = maxsize
        self.logger_name: str = WORKER_LOGGER_NAME
        self.router: Optional[logging.Handler] = None

    def setup(self, worker: Worker):
        """Setup dask worker plugin."""
        logger = logging.getLogger(self.logger_name)
        logger.setLevel(self.logging_level)
        logger.propagate = False # Root logger records handled by `logging.StreamHandler(sys.stdout)`
        self.router = MultiSocketHandler(logging_level=self.logging_level, maxsize=self.maxsize)
        self.router.setLevel(self.logging_level)
        # On dask workers, contextual information from `logging.LoggerAdapter` supersedes
        # these default filters (added for any root logger records from third-party dependencies)
        self.router.addFilter(DefaultProtocolNameFilter())
        self.router.addFilter(DefaultSocketAddressFilter())
        self.router.addFilter(DefaultTaskIdFilter())
        logger.addHandler(self.router)

    def teardown(self, worker: Worker):
        """Teardown dask worker plugin."""
        logger = logging.getLogger(self.logger_name)
        if self.router and self.router in logger.handlers:
            self.router.flush()
            logger.removeHandler(self.router)
            with suppress(Exception):
                self.router.close()
        self.router = None
