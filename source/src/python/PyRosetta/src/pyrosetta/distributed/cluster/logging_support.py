# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import logging
import os
import pickle
import select
import socketserver
import struct
import threading
import warnings

try:
    from distributed import WorkerPlugin
except ImportError:
    try:
        from distributed.diagnostics.plugin import WorkerPlugin
    except ImportError:
        print(
            "Importing 'pyrosetta.distributed.cluster.logging_support' requires the "
            + "third-party package 'distributed' as a dependency!\n"
            + "Please install the package into your python environment. "
            + "For installation instructions, visit:\n"
            + "https://pypi.org/project/distributed/\n"
        )
        raise

from contextlib import suppress
from functools import wraps
from typing import (
    Any,
    Callable,
    Generic,
    TypeVar,
    cast,
)


L = TypeVar("L", bound=Callable[..., Any])
G = TypeVar("G")


class LogRecordStreamHandler(socketserver.StreamRequestHandler):
    def handle(self):
        while True:
            chunk = self.connection.recv(4)
            if len(chunk) < 4:
                break
            slen = struct.unpack(">L", chunk)[0]
            chunk = self.connection.recv(slen)
            while len(chunk) < slen:
                chunk = chunk + self.connection.recv(slen - len(chunk))
            obj = pickle.loads(chunk)
            record = logging.makeLogRecord(obj)
            self.server.log_handler.handle(record)


class SocketListener(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    def __init__(
        self,
        host="localhost",
        port=logging.handlers.DEFAULT_TCP_LOGGING_PORT,
        handler=LogRecordStreamHandler,
    ):
        socketserver.ThreadingTCPServer.__init__(self, (host, port), LogRecordStreamHandler)
        self.abort = 0
        self.timeout = 1
        self.log_handler = handler
        self._thread = None

    def start(self):
        self._thread = threading.Thread(target=self.serve_forever, daemon=True)
        self._thread.start()

    def stop(self):
        self.shutdown()
        self.server_close()
        if self._thread:
            self._thread.join()

    def serve_until_stopped(self):
        abort = 0
        while not abort:
            rd, wr, ex = select.select([self.socket.fileno()], [], [], self.timeout)
            if rd:
                self.handle_request()
            abort = self.abort


class SocketLoggerPlugin(WorkerPlugin):
    def __init__(self, host, port):
        self.host = host
        self.port = port

    def setup(self, worker):
        worker.log_socket_address = (self.host, self.port)
        logger = logging.getLogger()
        logger.handlers.clear()
        logger.setLevel(logging.NOTSET)
        logger.addHandler(logging.handlers.SocketHandler(self.host, self.port))

    def teardown(self, worker):
        logging.shutdown()


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

        fh = logging.FileHandler(os.path.join(self.logs_path, "PyRosettaCluster.log",))
        fh.setFormatter(
            logging.Formatter(
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
        )
        logger.addHandler(fh)

    def _close_logger(self) -> None:
        """Close the logger for the client instance."""
        logger = logging.getLogger()
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            with suppress(Exception):
                handler.close()
        logging.shutdown()

    def _setup_socket_listeners(self, clients):
        for clients_index, client in clients.items():
            logging_file_split = list(os.path.splitext(self.logging_file))
            logging_file_split.insert(-1, f"_client-{clients_index}")
            logging_file = "".join(logging_file_split)
            host, port = client.run_on_scheduler(setup_socket_listener, logging_file, self.logging_level)
            plugin = SocketLoggerPlugin(host, port)
            plugin.idempotent = False # Always re-register plugin
            name = "socket-logger-plugin"
            if hasattr(client, "register_plugin"):
                client.register_plugin(plugin=plugin, name=name)
            else: # Deprecated in dask version 2023.9.2
                client.register_worker_plugin(plugin=plugin, name=name, nanny=False)

    def _close_socket_listeners(self, clients):
        for client in clients.values():
            client.run_on_scheduler(teardown_socket_listener)


def setup_socket_listener(logging_file, logging_level, dask_scheduler=None):
    logs_path = os.path.dirname(logging_file)
    if not os.path.isdir(logs_path):
        warnings.warn(
            f"Creating logs directory in 'LoggingSupport.setup_socket_listener': {logs_path}",
            RuntimeWarning,
            stacklevel=2,
        )
        os.makedirs(logs_path, exist_ok=True)

    handler = logging.FileHandler(logging_file, mode="a")
    handler.setLevel(logging_level)
    listener = SocketListener("0.0.0.0", 0, handler)
    listener.daemon = True
    listener.start()
    host, port = listener.socket.getsockname()
    dask_scheduler.log_listener = listener

    return host, port


def teardown_socket_listener(dask_scheduler=None):
    dask_scheduler.log_listener.stop()
    handler = dask_scheduler.log_listener.log_handler
    handler.flush()
    handler.close()
    logging.shutdown()


def setup_target_logging(func: L) -> L:
    """Support logging within the spawned thread."""
    @wraps(func)
    def wrapper(
        protocol,
        compressed_packed_pose,
        compressed_kwargs,
        q,
        logging_file,
        logging_level,
        log_socket_address,
        DATETIME_FORMAT,
        ignore_errors,
        protocols_key,
        decoy_ids,
        compression,
        client_residue_type_set,
        client_repr,
        **pyrosetta_init_kwargs,
    ):
        """Wrapper function to setup_target_logging."""
        logger = logging.getLogger()
        logger.setLevel(logging_level)
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            with suppress(Exception):
                handler.close()

        host, port = log_socket_address
        handler = logging.handlers.SocketHandler(host, port)
        logger.addHandler(handler)

        result = func(
            protocol,
            compressed_packed_pose,
            compressed_kwargs,
            q,
            logging_file,
            logging_level,
            log_socket_address,
            DATETIME_FORMAT,
            ignore_errors,
            protocols_key,
            decoy_ids,
            compression,
            client_residue_type_set,
            client_repr,
            **pyrosetta_init_kwargs,
        )

        logging.shutdown()

        return result

    return cast(L, wrapper)
