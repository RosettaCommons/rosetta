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
            self.server.handler.handle(record)


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
        self.handler = handler
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
    def __init__(self, host, port, logging_level):
        self.host = host
        self.port = port
        self.logging_level = logging_level

    def setup(self, worker):
        worker.socket_listener_address = (self.host, self.port)
        worker.logging_level = self.logging_level
        logger = logging.getLogger()
        logger.handlers.clear()
        logger.setLevel(self.logging_level)
        handler = logging.handlers.SocketHandler(self.host, self.port)
        logger.addHandler(handler)

    def teardown(self, worker):
        logger = logging.getLogger()
        for handler in logger.handlers[:]:
            handler.flush()
            logger.removeHandler(handler)
            with suppress(Exception):
                handler.close()


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

    def _setup_socket_listener(self, clients):
        logs_path = os.path.dirname(self.logging_file)
        if not os.path.isdir(logs_path):
            warnings.warn(
                f"Creating logs directory in 'LoggingSupport._setup_socket_listener': {logs_path}",
                RuntimeWarning,
                stacklevel=2,
            )
            os.makedirs(logs_path, exist_ok=True)

        handler = logging.FileHandler(self.logging_file, mode="a")
        handler.setLevel(self.logging_level)
        formatter = logging.Formatter(
            ":".join(
                [
                    "%(levelname)s",
                    # "%(protocol)s", # Extra key
                    "%(name)s",
                    "%(asctime)s",
                    " %(message)s",
                ]
            )
        )
        handler.setFormatter(formatter)
        self.socket_listener = SocketListener("0.0.0.0", 0, handler)
        self.socket_listener.daemon = True
        self.socket_listener.start()
        host, port = self.socket_listener.socket.getsockname()
        self._register_socket_logger_plugins(clients, host, port)

    def _register_socket_logger_plugins(self, clients, host, port):
        for client in clients.values():
            plugin = SocketLoggerPlugin(host, port, self.logging_level)
            plugin.idempotent = False # Always re-register plugin
            name = "socket-logger-plugin"
            if hasattr(client, "register_plugin"):
                client.register_plugin(plugin=plugin, name=name)
            else: # Deprecated since dask version 2023.9.2
                client.register_worker_plugin(plugin=plugin, name=name, nanny=False)

    def _close_socket_listener(self):
        self.socket_listener.stop()
        handler = self.socket_listener.handler
        handler.flush()
        with suppress(Exception):
            handler.close()

def setup_target_logging(func: L) -> L:
    """Support logging within the spawned thread."""
    @wraps(func)
    def wrapper(
        protocol,
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
        **pyrosetta_init_kwargs,
    ):
        """Wrapper function to setup_target_logging."""
        logger = logging.getLogger()
        logger.setLevel(logging_level)
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            with suppress(Exception):
                handler.close()

        host, port = socket_listener_address
        handler = logging.handlers.SocketHandler(host, port)
        logger.addHandler(handler)

        result = func(
            protocol,
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
            **pyrosetta_init_kwargs,
        )

        handler.flush()
        logger.removeHandler(handler)
        with suppress(Exception):
            handler.close()

        return result

    return cast(L, wrapper)
