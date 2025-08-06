# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    import billiard
    import msgpack
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.logging_support' requires the "
        + "third-party packages 'billiard' and 'msgpack' as dependencies!\n"
        + "Please install the packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/billiard/\n"
        + "https://pypi.org/project/msgpack/\n"
    )
    raise

import hashlib
import hmac
import logging
import os
import select
import socketserver
import struct
import threading
import traceback
import warnings

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
            header = self.connection.recv(4)
            if len(header) < 4:
                break
            msg_len = struct.unpack(">L", header)[0]
            if msg_len <= 0 or msg_len > self.server.max_packet_size:
                break
            msg = self.connection.recv(msg_len)
            while len(msg) < msg_len:
                msg += self.connection.recv(msg_len - len(msg))
            try:
                obj = self.unPickle(msg)
                record = logging.makeLogRecord(obj)
                self.server.handler.handle(record)
            except Exception as ex:
                warnings.warn(
                    f"{type(ex).__name__}: {ex}. Rejected log packet:\n{traceback.format_exc()}",
                    RuntimeWarning,
                    stacklevel=2,
                )

    def unPickle(self, msg: bytes) -> Dict[str, Any]:
        packet = msgpack.unpackb(msg, raw=False)
        signature = packet["signature"]
        compressed_record = packet["compressed_record"]
        required_signature = hmac.new(self.server.hmac_key, compressed_record, hashlib.sha256).digest()
        if not hmac.compare_digest(required_signature, signature):
            raise ValueError("Logging socket listener received a bad hash-based message authentication code!")

        return msgpack.unpackb(compressed_record, raw=False)


class SocketListener(socketserver.ThreadingTCPServer):
    """
    TCP socket-based logging receiver modified from logging cookbook recipe:
    https://docs.python.org/3/howto/logging-cookbook.html#sending-and-receiving-logging-events-across-a-network
    """
    allow_reuse_address = True
    def __init__(self, host: str, port: int, handler: logging.Handler, hmac_key: bytes, timeout: Union[float, int]) -> None:
        super().__init__((host, port), LogRecordRequestHandler)
        self.handler = handler
        self.hmac_key = hmac_key
        self.max_packet_size = 10 * 1024 * 1024 # Maximum of 10 MiB per log message
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


class MaskedBytes(bytes):
    """A `bytes` subclass to mask its contents if the `repr` method is called."""
    def __new__(cls, value: bytes) -> bytes:
        return super().__new__(cls, value)

    def __repr__(self) -> str:
        return "<masked-value>"


class MsgpackHmacSocketHandler(logging.handlers.SocketHandler):
    _supported_type = (str, int, float, bool, type(None))

    def __init__(self, host: str, port: int, hmac_key: bytes) -> None:
        super().__init__(host, port)
        self.hmac_key = hmac_key

    def makePickle(self, record: logging.LogRecord) -> bytes:
        """Compress a logging record with MessagePack and hash-based message authentication codes (HMAC)."""
        record_dict = dict(
            msg=record.msg,
            args=record.args,
            name=record.name,
            levelno=int(record.levelno),
            levelname=record.levelname,
            pathname=record.pathname,
            lineno=int(record.lineno),
            module=record.module,
            funcName=record.funcName,
            created=float(record.created),
            msecs=float(record.msecs),
            relativeCreated=float(record.relativeCreated),
            process=int(record.process),
            processName=record.processName,
            thread=int(record.thread),
            threadName=record.threadName,
            protocol_name=record.protocol_name,
        )
        if record.exc_info:
            record_dict["exc_text"] = logging.Formatter().formatException(record.exc_info)
        if record.stack_info:
            record_dict["stack_info"] = record.stack_info

        compressed_record = msgpack.packb(record_dict, use_bin_type=True)
        signature = hmac.new(self.hmac_key, compressed_record, hashlib.sha256).digest()
        packet = dict(signature=signature, compressed_record=compressed_record, version=1.0)
        compressed_packet = msgpack.packb(packet, use_bin_type=True)
        header = struct.pack(">L", len(compressed_packet))

        return header + compressed_packet


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
        masked_key = MaskedBytes(os.urandom(32))
        _host, _port = tuple(s.strip() for s in self.logging_address.split(":"))
        self.socket_listener = SocketListener(_host, int(_port), handler, masked_key, self.timeout)
        self.socket_listener.daemon = True
        self.socket_listener.start()
        socket_listener_address = self.socket_listener.socket.getsockname()
        host, port = socket_listener_address
        logging.info(f"Logging socket listener: http://{host}:{port}")

        return socket_listener_address, masked_key

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
        # User-provided PyRosetta protocol name is the first argument of 'user_spawn_thread' and 'target'
        protocol_name = args[0]
        value = current_protocol_name.set(protocol_name)
        try:
            return func(*args, **kwargs)
        finally:
            current_protocol_name.reset(value)

    return cast(L, wrapper)


def setup_logger(
    logging_level: str, socket_listener_address: Tuple[str, int], masked_key: bytes
) -> Tuple[logging.RootLogger, logging.handlers.SocketHandler]:
    """Setup socket logging handler."""
    logger = logging.getLogger()
    logger.setLevel(logging_level)
    for _handler in logger.handlers[:]:
        for _filter in _handler.filters[:]:
            if isinstance(_filter, ProtocolContextFilter):
                _handler.removeFilter(_filter)
        logger.removeHandler(_handler)
        with suppress(Exception):
            _handler.close()

    host, port = socket_listener_address
    socket_handler = MsgpackHmacSocketHandler(host, port, masked_key)
    context_filter = ProtocolContextFilter()
    socket_handler.addFilter(context_filter)
    socket_handler.closeOnError = True
    logger.addHandler(socket_handler)

    return logger, socket_handler, context_filter


def close_logger(
    logger: logging.RootLogger,
    socket_handler: logging.handlers.SocketHandler,
    context_filter: ProtocolContextFilter,
) -> None:
    """Teardown socket logging handler."""
    socket_handler.flush()
    socket_handler.removeFilter(context_filter)
    logger.removeHandler(socket_handler)
    with suppress(Exception):
        socket_handler.close()


def setup_target_logging(func: L) -> L:
    """Support logging within the billiard spawned thread."""
    @bind_protocol
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
        logger, socket_handler, context_filter = setup_logger(
            logging_level, socket_listener_address, masked_key
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

        close_logger(logger, socket_handler, context_filter)

        return result

    return cast(L, wrapper)


def setup_worker_logging(func: L) -> L:
    """Support logging within the dask worker thread."""
    @bind_protocol
    @wraps(func)
    def wrapper(
        protocol_name: str,
        compressed_protocol: bytes,
        compressed_packed_pose: bytes,
        compressed_kwargs: bytes,
        pyrosetta_init_kwargs: Dict[str, Any],
        client_repr: str,
        extra_args: Dict[str, Any],
    ):
        """Wrapper function to setup_worker_logging."""
        logging_level = extra_args["logging_level"]
        socket_listener_address = extra_args["socket_listener_address"]
        masked_key = extra_args["masked_key"]

        logger, socket_handler, context_filter = setup_logger(
            logging_level, socket_listener_address, masked_key
        )

        result = func(
            protocol_name,
            compressed_protocol,
            compressed_packed_pose,
            compressed_kwargs,
            pyrosetta_init_kwargs,
            client_repr,
            extra_args,
        )

        close_logger(logger, socket_handler, context_filter)

        return result

    return cast(L, wrapper)
