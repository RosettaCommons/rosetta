# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    import billiard
    from distributed import Client, Worker
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.logging_support' requires the "
        + "third-party packages 'billiard' and 'distributed' as dependencies!\n"
        + "Please install the packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/billiard/\n"
        + "https://pypi.org/project/distributed/\n"
    )
    raise

import logging
import os
import warnings

from contextlib import suppress
from functools import wraps
from typing import (
    AbstractSet,
    Any,
    Callable,
    Dict,
    Generic,
    List,
    Optional,
    Tuple,
    TypeVar,
    Union,
    cast,
)

from pyrosetta.distributed.cluster.worker_plugins import (
    SocketLoggerPlugin,
    SOCKET_LOGGER_PLUGIN_NAME,
    WORKER_LOGGER_NAME,
)
from pyrosetta.distributed.cluster.logging_filters import (
    SetProtocolNameFilter,
    SetSocketAddressFilter,
    format_socket_address,
    split_socket_address,
)
from pyrosetta.distributed.cluster.logging_handlers import MsgpackHmacSocketHandler
from pyrosetta.distributed.cluster.logging_listeners import SocketListener


G = TypeVar("G")
L = TypeVar("L", bound=Callable[..., Any])
Q = TypeVar("Q", bound=billiard.Queue)


class LoggingSupport(Generic[G]):
    """Supporting logging methods for PyRosettaCluster."""
    def __init__(self) -> None:
        """Log warnings from the warnings module."""
        logging.captureWarnings(True)

    def _setup_logger(self) -> None:
        """Open the logger for the client instance."""
        logger = logging.getLogger()
        logger.setLevel(self.logging_level)
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            with suppress(Exception):
                handler.close()

        if not os.path.isdir(self.logs_path):
            warnings.warn(
                f"Creating logs directory in 'LoggingSupport._setup_logger': {self.logs_path}",
                RuntimeWarning,
                stacklevel=2,
            )
            os.makedirs(self.logs_path, exist_ok=True)

        handler = logging.FileHandler(os.path.join(self.logs_path, "PyRosettaCluster.log"))
        formatter = logging.Formatter(
            ":".join(
                [
                    "PyRosettaCluster",
                    "%(asctime)s",
                    "%(levelname)s",
                    "%(name)s",
                    " %(message)s",
                ]
            )
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    def _close_logger(self) -> None:
        """Close the logger for the client instance."""
        logger = logging.getLogger()
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            with suppress(Exception):
                handler.close()
        logging.shutdown()

    def _setup_socket_listener(self, clients: Dict[int, Client]) -> Tuple[Tuple[str, int], bytes]:
        """Setup logging socket listener."""
        logs_path = os.path.dirname(self.logging_file)
        if not os.path.isdir(logs_path):
            warnings.warn(
                f"Creating logs directory in 'LoggingSupport._setup_socket_listener': {logs_path}",
                RuntimeWarning,
                stacklevel=2,
            )
            os.makedirs(logs_path, exist_ok=True)

        self.socket_listener = SocketListener(
            self.logging_address,
            self.logging_file,
            self.logging_level,
            self.timeout,
            self.ignore_errors,
        )
        self.socket_listener.daemon = True
        self.socket_listener.start()
        socket_listener_address = self.socket_listener.socket_listener_address
        masked_key = self.socket_listener.masked_key
        logging.info("Logging socket listener: http://{0}:{1}".format(*socket_listener_address))
        self._register_socket_logger_plugin(clients)

        return socket_listener_address, masked_key

    def _register_socket_logger_plugin(self, clients: Dict[int, Client]) -> None:
        """Register `SocketLoggerPlugin` as a dask worker plugin on dask clients."""
        for client in clients.values():
            plugin = SocketLoggerPlugin(self.logging_level, maxsize=32)
            plugin.idempotent = True # Never re-register plugin
            if hasattr(client, "register_plugin"):
                client.register_plugin(plugin=plugin, name=SOCKET_LOGGER_PLUGIN_NAME)
            else: # Deprecated since dask version 2023.9.2
                client.register_worker_plugin(plugin=plugin, name=SOCKET_LOGGER_PLUGIN_NAME, nanny=False)

    def _close_socket_listener(self, clients: Dict[int, Client]) -> None:
        """Close logging socket listener."""
        self._close_socket_logger_plugins(clients)
        self.socket_listener.stop()
        handler = self.socket_listener.handler
        handler.flush()
        with suppress(Exception):
            handler.close()

    def _close_socket_logger_plugins(self, clients: Dict[int, Client]) -> None:
        """Purge cached logging socket addresses on all dask workers."""
        socket_listener_address = self.socket_listener.socket_listener_address
        for client in clients.values():
            results = client.run(
                self._purge_socket_logger_plugin_address,
                socket_listener_address,
                workers=None,
                wait=True,
                nanny=False,
                on_error="return",
            )
            for worker_address, result in results.items():
                if result:
                    logging.warning(
                        f"Logger was not closed cleanly on dask worker ({worker_address}) - {result}"
                    )

    def _purge_socket_logger_plugin_address(
        self,
        socket_listener_address: Tuple[str, int],
        dask_worker: Worker,
    ) -> None:
        """Close and remove an item from the worker logger plugin router."""
        router = dask_worker.plugins[SOCKET_LOGGER_PLUGIN_NAME]._router
        with suppress(Exception):
            router.purge_address(socket_listener_address)


def setup_target_logger(
    protocol_name: str,
    socket_listener_address: Tuple[str, int],
    masked_key: bytes,
    logging_level: str,
) -> Tuple[logging.RootLogger, logging.handlers.SocketHandler, List[logging.Filter]]:
    """Setup socket logging handler."""
    logger = logging.getLogger()
    logger.setLevel(logging_level)
    for _handler in logger.handlers[:]:
        for _filter in _handler.filters[:]:
            if isinstance(_filter, (SetProtocolNameFilter, SetSocketAddressFilter)):
                _handler.removeFilter(_filter)
        logger.removeHandler(_handler)
        with suppress(Exception):
            _handler.close()

    host, port = socket_listener_address
    socket_handler = MsgpackHmacSocketHandler(host, port, masked_key=masked_key)
    filters = [
        SetProtocolNameFilter(protocol_name),
        SetSocketAddressFilter(socket_listener_address),
    ]
    for _filter in filters:
        socket_handler.addFilter(_filter)
    socket_handler.closeOnError = True
    logger.addHandler(socket_handler)

    return logger, socket_handler, filters


def close_target_logger(
    logger: logging.RootLogger,
    socket_handler: logging.handlers.SocketHandler,
    filters: List[logging.Filter],
) -> None:
    """Teardown socket logging handler."""
    socket_handler.flush()
    for _filter in filters:
        socket_handler.removeFilter(_filter)
    logger.removeHandler(socket_handler)
    with suppress(Exception):
        socket_handler.close()


def setup_target_logging(func: L) -> L:
    """Support logging within the billiard spawned thread."""
    @wraps(func)
    def wrapper(
        protocol_name: str,
        compressed_protocol: bytes,
        compressed_packed_pose: bytes,
        compressed_kwargs: bytes,
        q: Q,
        logging_level: str,
        socket_listener_address: Tuple[str, int],
        datetime_format: str,
        ignore_errors: bool,
        protocols_key: str,
        decoy_ids: List[int],
        compression: Optional[Union[str, bool]],
        client_residue_type_set: AbstractSet[str],
        client_repr: str,
        masked_key: bytes,
        **pyrosetta_init_kwargs: Dict[str, Any],
    ):
        """Wrapper function to setup_target_logging."""
        logger, socket_handler, filters = setup_target_logger(
            protocol_name, socket_listener_address, masked_key, logging_level
        )

        result = func(
            protocol_name,
            compressed_protocol,
            compressed_packed_pose,
            compressed_kwargs,
            q,
            logging_level,
            socket_listener_address,
            datetime_format,
            ignore_errors,
            protocols_key,
            decoy_ids,
            compression,
            client_residue_type_set,
            client_repr,
            masked_key,
            **pyrosetta_init_kwargs,
        )

        close_target_logger(logger, socket_handler, filters)

        return result

    return cast(L, wrapper)
