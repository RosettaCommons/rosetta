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

from pyrosetta.distributed.cluster.logging_handlers import MultiSocketHandler
from pyrosetta.distributed.cluster.logging_filters import (
    DefaultProtocolNameFilter,
    DefaultSocketAddressFilter,
)


SOCKET_LOGGER_PLUGIN_NAME: str = "PyRosettaCluster_socket_logger_plugin"
WORKER_LOGGER_NAME: str = "PyRosettaCluster_dask_worker"


class SocketLoggerPlugin(WorkerPlugin):
    """Install a `MultiSocketHandler` logging handler on a dask worker logger."""
    def __init__(
        self,
        logging_level: str,
        maxsize=128,
    ):
        self.logging_level = logging_level
        self.maxsize = maxsize
        self.logger_name = WORKER_LOGGER_NAME
        self._router = None

    def setup(self, worker: Worker):
        """Setup dask worker plugin."""
        logger = logging.getLogger(self.logger_name)
        logger.setLevel(self.logging_level)
        self._router = MultiSocketHandler(maxsize=self.maxsize)
        # On dask workers, contextual information from `logging.LoggerAdapter` supersedes
        # these default filters (added for any third-party software root loggers)
        self._router.addFilter(DefaultProtocolNameFilter())
        self._router.addFilter(DefaultSocketAddressFilter())
        logger.addHandler(self._router)

    def teardown(self, worker: Worker):
        """Teardown dask worker plugin."""
        logger = logging.getLogger(self.logger_name)
        if self._router and self._router in logger.handlers:
            self._router.flush()
            logger.removeHandler(self._router)
            with suppress(Exception):
                self._router.close()
        self._router = None
