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
    from distributed import Client
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.logging_support' requires the "
        + "third-party packages 'billiard', 'msgpack', and 'distributed' as dependencies!\n"
        + "Please install the packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/billiard/\n"
        + "https://pypi.org/project/msgpack/\n"
        + "https://pypi.org/project/distributed/\n"
    )
    raise

import collections
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
from functools import wraps
from typing import (
    AbstractSet,
    Any,
    Callable,
    Dict,
    Generic,
    List,
    Optional,
    OrderedDict,
    Tuple,
    TypeVar,
    Union,
    cast,
)


G = TypeVar("G")
L = TypeVar("L", bound=Callable[..., Any])
Q = TypeVar("Q", bound=billiard.Queue)


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
                _err_msg = f"{type(ex).__name__}: {ex}. Rejected log packet:\n{traceback.format_exc()}"
                if self.server.ignore_errors:
                    warnings.warn(_err_msg, RuntimeWarning, stacklevel=2)
                else:
                    raise BufferError(_err_msg)

    def unPickle(self, msg: bytes) -> Dict[str, Any]:
        packet = msgpack.unpackb(msg, raw=False)
        signature = packet["signature"]
        compressed_record = packet["compressed_record"]
        required_signature = hmac.new(self.server.masked_key, compressed_record, hashlib.sha256).digest()
        if not hmac.compare_digest(required_signature, signature):
            raise ValueError("Logging socket listener received a bad hash-based message authentication code!")

        return msgpack.unpackb(compressed_record, raw=False)


class SocketListener(socketserver.ThreadingTCPServer):
    """
    TCP socket-based logging receiver modified from logging cookbook recipe:
    https://docs.python.org/3/howto/logging-cookbook.html#sending-and-receiving-logging-events-across-a-network
    """
    allow_reuse_address = True

    def __init__(
        self,
        logging_address: str,
        logging_file: str,
        logging_level: str,
        timeout: Union[float, int],
        ignore_errors: bool,
    ) -> None:
        _host, _port = tuple(s.strip() for s in logging_address.split(":"))
        super().__init__((_host, int(_port)), LogRecordRequestHandler)
        self.socket_listener_address = self.socket.getsockname()
        self.handler = self.setup_handler(logging_file, logging_level)
        self.masked_key = MaskedBytes(os.urandom(32))
        self.timeout = timeout
        self.ignore_errors = ignore_errors
        self.max_packet_size = 10 * 1024 * 1024 # Maximum of 10 MiB per log message
        self.abort = 0
        self._thread = None

    def setup_handler(self, logging_file: str, logging_level: str) -> logging.FileHandler:
        """Setup logging file handler for logging socket listener."""
        handler = logging.FileHandler(logging_file, mode="a")
        handler.setLevel(logging_level)
        handler.addFilter(SocketAddressFilter(self.socket_listener_address))
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

        return handler

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
        return "mask"


class MsgpackHmacSocketHandler(logging.handlers.SocketHandler):
    _supported_type = (str, int, float, bool, type(None))

    def __init__(self, host: str, port: int, masked_key: bytes) -> None:
        super().__init__(host, port)
        self.masked_key = masked_key

    def makePickle(self, record: logging.LogRecord) -> bytes:
        """Compress a logging record with MessagePack and a hash-based message authentication code (HMAC)."""
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
            protocol_name=record.protocol_name, # Set by `SetProtocolNameFilter``
            socket_address=record.socket_address, # Set by `SetSocketAddressFilter`
        )
        if record.exc_info:
            record_dict["exc_text"] = logging.Formatter().formatException(record.exc_info)
        if record.stack_info:
            record_dict["stack_info"] = record.stack_info

        compressed_record = msgpack.packb(record_dict, use_bin_type=True)
        signature = hmac.new(self.masked_key, compressed_record, hashlib.sha256).digest()
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

        return socket_listener_address, masked_key

    def _close_socket_listener(self) -> None:
        """Close logging socket listener."""
        self.socket_listener.stop()
        handler = self.socket_listener.handler
        handler.flush()
        with suppress(Exception):
            handler.close()

    def _close_worker_loggers(self, clients: Dict[int, Client]) -> None:
        """Closer dask worker loggers."""
        for client in clients.values():
            _ = client.run(
                close_worker_loggers,
                workers=None,
                wait=True,
                nanny=False,
                on_error="ignore",
            )


class SetProtocolNameFilter(logging.Filter):
    """Set protocol name for logging socket listener formatter."""
    def __init__(self, protocol_name: str) -> None:
        super().__init__()
        self.protocol_name = protocol_name

    def filter(self, record: logging.LogRecord) -> bool:
        record.protocol_name = self.protocol_name
        return True


class SetSocketAddressFilter(logging.Filter):
    """Set socket address for logging socket listener filter."""
    def __init__(self, socket_listener_address: Tuple[str, int]) -> None:
        super().__init__()
        self.socket_address = _norm_socket_address(socket_listener_address)

    def filter(self, record: logging.LogRecord) -> bool:
        record.socket_address = self.socket_address
        return True


class SocketAddressFilter(logging.Filter):
    """Filter log records for the logging socket listener address."""
    def __init__(self, socket_listener_address: Tuple[str, int]) -> None:
        super().__init__()
        self.socket_address = _norm_socket_address(socket_listener_address)

    def filter(self, record: logging.LogRecord) -> bool:
        return record.socket_address == self.socket_address


class WorkerLoggerLRUCache(Generic[G]):
    """
    Memoize dask worker loggers up to a maximum sized cache, pruning least recently used (LRU)
    dask worker loggers first.
    """
    def __init__(self, maxsize: int = 128) -> None:
        self.cache: OrderedDict[
            Tuple[str, Tuple[str, int]],
            Tuple[logging.RootLogger, logging.handlers.SocketHandler, List[logging.Filter]]
        ] = collections.OrderedDict()
        self.maxsize: int = maxsize

    def to_key(
        self, protocol_name: str, socket_listener_address: Tuple[str, int]
    ) -> Tuple[str, Tuple[str, int]]:
        """Get normalized key for worker logger cache."""
        return (protocol_name, socket_listener_address)

    def maybe_prune(self) -> None:
        """Prune worker logger cache to maximum capacity."""
        if len(self.cache) >= self.maxsize:
            _, (logger, socket_handler, filters) = self.cache.popitem(last=False)
            close_logger(logger, socket_handler, filters)

    def put(
        self,
        protocol_name: str,
        socket_listener_address: Tuple[str, int],
        masked_key: bytes,
        logging_level: str,
    ) -> None:
        """Add item to worker logger cache if not already added, and prune worker logger cache."""
        key = self.to_key(protocol_name, socket_listener_address)
        if key not in self.cache:
            self.cache[key] = setup_logger(
                protocol_name, socket_listener_address, masked_key, logging_level
            )
        self.cache.move_to_end(key, last=True)
        self.maybe_prune()

    def clear(self) -> None:
        """Close worker loggers and clear worker logger cache."""
        for logger, socket_handler, filters in self.cache.values():
            close_logger(logger, socket_handler, filters)
        self.cache.clear()


# Instantiate worker logger cache in module scope for worker imports
worker_logger_cache: WorkerLoggerLRUCache = WorkerLoggerLRUCache(maxsize=32)


def _norm_socket_address(socket_listener_address: Tuple[str, int]) -> str:
    """Normalize a socket listener address for socket listener handler filters."""
    return ":".join(map(str, socket_listener_address))


def setup_logger(
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
    socket_handler = MsgpackHmacSocketHandler(host, port, masked_key)
    filters = [
        SetProtocolNameFilter(protocol_name),
        SetSocketAddressFilter(socket_listener_address),
    ]
    for _filter in filters:
        socket_handler.addFilter(_filter)
    socket_handler.closeOnError = True
    logger.addHandler(socket_handler)

    return logger, socket_handler, filters


def close_logger(
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
        logger, socket_handler, filters = setup_logger(
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

        close_logger(logger, socket_handler, filters)

        return result

    return cast(L, wrapper)


def setup_worker_logging(func: L) -> L:
    """Support logging within the dask worker thread."""
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
        worker_logger_cache.put(protocol_name, socket_listener_address, masked_key, logging_level)

        return func(
            protocol_name,
            compressed_protocol,
            compressed_packed_pose,
            compressed_kwargs,
            pyrosetta_init_kwargs,
            client_repr,
            extra_args,
        )

    return cast(L, wrapper)


def close_worker_loggers() -> None:
   """Clear dask worker logger cache."""
   worker_logger_cache.clear()
