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
    from distributed import Client, Worker, WorkerPlugin
except ImportError:
    try:
        from distributed import Client, Worker
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
from contextvars import ContextVar
from functools import wraps
from typing import (
    Any,
    Callable,
    Dict,
    Generic,
    TypeVar,
    Union,
    cast,
)


L = TypeVar("L", bound=Callable[..., Any])
G = TypeVar("G")

SOCKET_LOGGER_PLUGIN_NAME: str = "socket-logger-plugin"
DEFAULT_PROTOCOL_NAME: str = "user_spawn_thread"


class LogRecordRequestHandler(socketserver.StreamRequestHandler):
    """
    Handler for a streaming logging request modified from logging cookbook recipe:
    https://docs.python.org/3/howto/logging-cookbook.html#sending-and-receiving-logging-events-across-a-network
    """
    def handle(self) -> None:
        while True:
            chunk = self.connection.recv(4)
            if len(chunk) < 4:
                break
            slen = struct.unpack(">L", chunk)[0]
            chunk = self.connection.recv(slen)
            while len(chunk) < slen:
                chunk = chunk + self.connection.recv(slen - len(chunk))
            obj = self.decompress(chunk)
            record = logging.makeLogRecord(obj)
            self.server.handler.handle(record)

    def decompress(self, chunk: bytes) -> Any:
        return pickle.loads(chunk)


class SocketListener(socketserver.ThreadingTCPServer):
    """
    TCP socket-based logging receiver modified from logging cookbook recipe:
    https://docs.python.org/3/howto/logging-cookbook.html#sending-and-receiving-logging-events-across-a-network
    """
    allow_reuse_address = True
    def __init__(self, host: str, port: int, handler: logging.Handler, timeout: Union[float, int]) -> None:
        super().__init__((host, port), LogRecordRequestHandler)
        self.handler = handler
        self.timeout = timeout
        self.abort = 0
        self._thread = None

    def start(self) -> None:
        self._thread = threading.Thread(target=self.serve_forever, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self.shutdown()
        self.server_close()
        if self._thread:
            self._thread.join()

    def serve_until_stopped(self) -> None:
        abort = 0
        while not abort:
            _read_ready, _, _ = select.select([self.socket.fileno()], [], [], self.timeout)
            if _read_ready:
                self.handle_request()
            abort = self.abort


class SocketLoggerPlugin(WorkerPlugin):
    """Dask worker plugin for logging socket handler."""
    def __init__(self, host: str, port: int, logging_level: str) -> None:
        self.host = host
        self.port = port
        self.logging_level = logging_level

    def setup(self, worker: Worker) -> None:
        """Setup dask worker plugin for logging socket handler."""
        logger = logging.getLogger()
        logger.handlers.clear()
        logger.setLevel(self.logging_level)
        handler = logging.handlers.SocketHandler(self.host, self.port)
        handler.addFilter(ProtocolContextFilter())
        handler.closeOnError = True
        logger.addHandler(handler)

    def teardown(self, worker: Worker) -> None:
        """Teardown dask worker plugin for logging socket handler."""
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

    def _setup_socket_listener(self, clients: Dict[int, Client]) -> None:
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
                    "%(protocol)s", # Extra key
                    "%(name)s",
                    "%(asctime)s",
                    " %(message)s",
                ]
            )
        )
        handler.setFormatter(formatter)
        handler.addFilter(ProtocolDefaultFilter())
        _host, _port = map(lambda s: s.strip(), self.logging_address.split(":"))
        self.socket_listener = SocketListener(_host, int(_port), handler, self.timeout)
        self.socket_listener.daemon = True
        self.socket_listener.start()
        host, port = self.socket_listener.socket.getsockname()
        logging.info(f"Logging socket listener: http://{host}:{port}")
        self._register_socket_logger_plugins(clients, host, port)

    def _register_socket_logger_plugins(self, clients: Dict[int, Client], host: str, port: int) -> None:
        for client in clients.values():
            plugin = SocketLoggerPlugin(host, port, self.logging_level)
            plugin.idempotent = False # Always re-register plugin
            if hasattr(client, "register_plugin"):
                client.register_plugin(plugin=plugin, name=SOCKET_LOGGER_PLUGIN_NAME)
            else: # Deprecated since dask version 2023.9.2
                client.register_worker_plugin(plugin=plugin, name=SOCKET_LOGGER_PLUGIN_NAME, nanny=False)

    def _close_socket_listener(self) -> None:
        self.socket_listener.stop()
        handler = self.socket_listener.handler
        handler.flush()
        with suppress(Exception):
            handler.close()


class ProtocolDefaultFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        if not hasattr(record, "protocol"):
            record.protocol = DEFAULT_PROTOCOL_NAME
        return True


class ProtocolContextFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        record.protocol = current_protocol.get()
        return True


current_protocol = ContextVar("current_protocol", default=DEFAULT_PROTOCOL_NAME)


def bind_protocol(func: L) -> L:
    @wraps(func)
    def wrapper(*args, **kwargs):
        # User-provided PyRosetta protocol is the first argument of 'user_spawn_thread' and 'target'
        protocol = args[0]
        value = current_protocol.set(protocol.__name__)
        try:
            return func(*args, **kwargs)
        finally:
            current_protocol.reset(value)

    return cast(L, wrapper)


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
        handler.addFilter(ProtocolContextFilter())
        handler.closeOnError = True
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
