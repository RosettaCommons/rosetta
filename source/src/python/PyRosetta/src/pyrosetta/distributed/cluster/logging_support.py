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
    import billiard
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.logging_support' requires the "
        + "third-party package 'billiard' as a dependency!\n"
        + "Please install the package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/billiard/\n"
    )
    raise

from contextlib import suppress
from contextvars import ContextVar
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


G = TypeVar("G")
L = TypeVar("L", bound=Callable[..., Any])
Q = TypeVar("Q", bound=billiard.Queue)


DEFAULT_PROTOCOL_NAME: str = "PyRosettaCluster"


class LogRecordRequestHandler(socketserver.StreamRequestHandler):
    """
    Handler for a streaming logging request modified from logging cookbook recipe:
    https://docs.python.org/3/howto/logging-cookbook.html#sending-and-receiving-logging-events-across-a-network
    """
    def handle(self) -> None:
        while True:
            msg = self.connection.recv(4)
            if len(msg) < 4:
                break
            msg_len = struct.unpack(">L", msg)[0]
            msg = self.connection.recv(msg_len)
            while len(msg) < msg_len:
                msg += self.connection.recv(msg_len - len(msg))
            obj = self.decompress(msg)
            record = logging.makeLogRecord(obj)
            self.server.handler.handle(record)

    def decompress(self, msg: bytes) -> Any:
        return pickle.loads(msg)


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

    def _setup_socket_listener(self) -> Tuple[str, int]:
        """Setup socket logging socket listener."""
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
                    "%(protocol_name)s", # Extra key
                    "%(name)s",
                    "%(asctime)s",
                    " %(message)s",
                ]
            )
        )
        handler.setFormatter(formatter)
        handler.addFilter(ProtocolDefaultFilter())
        _host, _port = iter(s.strip() for s in self.logging_address.split(":"))
        self.socket_listener = SocketListener(_host, int(_port), handler, self.timeout)
        self.socket_listener.daemon = True
        self.socket_listener.start()
        host, port = self.socket_listener.socket.getsockname()
        logging.info(f"Logging socket listener: http://{host}:{port}")

        return host, port

    def _close_socket_listener(self) -> None:
        """Teardown socket logging socket listener."""
        self.socket_listener.stop()
        handler = self.socket_listener.handler
        handler.flush()
        with suppress(Exception):
            handler.close()


class ProtocolDefaultFilter(logging.Filter):
    """Set default protocol name for logging socket listener formatter."""
    def filter(self, record: logging.LogRecord) -> bool:
        if not hasattr(record, "protocol_name"):
            record.protocol_name = DEFAULT_PROTOCOL_NAME
        return True


class ProtocolContextFilter(logging.Filter):
    """Set bound protocol name for logging socket listener formatter."""
    def filter(self, record: logging.LogRecord) -> bool:
        record.protocol_name = current_protocol_name.get()
        return True


current_protocol_name: ContextVar = ContextVar("current_protocol_name", default=DEFAULT_PROTOCOL_NAME)


def bind_protocol(func: L) -> L:
    """Set and reset current protocol name for socket logging filters."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        """Wrapper function to bind_protocol."""
        # User-provided PyRosetta protocol is the first argument of 'user_spawn_thread' and 'target'
        protocol = args[0]
        value = current_protocol_name.set(protocol.__name__)
        try:
            return func(*args, **kwargs)
        finally:
            current_protocol_name.reset(value)

    return cast(L, wrapper)


def setup_logger(
    logging_level: str, socket_listener_address: Tuple[str, int]
) -> Tuple[logging.RootLogger, logging.handlers.SocketHandler]:
    """Setup socket logging handler."""
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

    return logger, handler


def close_logger(
    logger: logging.RootLogger, handler: logging.handlers.SocketHandler
) -> None:
    """Teardown socket logging handler."""
    handler.flush()
    logger.removeHandler(handler)
    with suppress(Exception):
        handler.close()


def setup_target_logging(func: L) -> L:
    """Support logging within the billiard spawned thread."""
    @bind_protocol
    @wraps(func)
    def wrapper(
        protocol: Callable[..., Any],
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
        **pyrosetta_init_kwargs: Dict[str, Any],
    ):
        """Wrapper function to setup_target_logging."""
        logger, handler = setup_logger(logging_level, socket_listener_address)

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

        close_logger(logger, handler)

        return result

    return cast(L, wrapper)


def setup_worker_logging(func: L) -> L:
    """Support logging within the dask worker thread."""
    @bind_protocol
    @wraps(func)
    def wrapper(
        protocol: Callable[..., Any],
        compressed_packed_pose: bytes,
        compressed_kwargs: bytes,
        pyrosetta_init_kwargs: Dict[str, Any],
        extra_args: Dict[str, Any],
    ):
        """Wrapper function to setup_worker_logging."""
        logging_level = extra_args["logging_level"]
        socket_listener_address = extra_args["socket_listener_address"]

        logger, handler = setup_logger(logging_level, socket_listener_address)

        result = func(
            protocol,
            compressed_packed_pose,
            compressed_kwargs,
            pyrosetta_init_kwargs,
            extra_args,
        )

        close_logger(logger, handler)

        return result

    return cast(L, wrapper)
