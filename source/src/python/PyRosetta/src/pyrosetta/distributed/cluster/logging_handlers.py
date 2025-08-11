# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    import msgpack
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.logging_handlers' requires the "
        + "third-party package 'msgpack' as a dependency!\n"
        + "Please install the package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/msgpack/\n"
    )
    raise

import collections
import logging
import struct
import sys
import types

from contextlib import contextmanager, suppress
from functools import wraps
from typing import (
    Any,
    Callable,
    Dict,
    Generator,
    Optional,
    OrderedDict,
    Tuple,
    TypeVar,
    cast,
)

from pyrosetta.distributed.cluster.hkdf import hmac_digest
from pyrosetta.distributed.cluster.logging_filters import split_socket_address


L = TypeVar("L", bound=Callable[..., Any])


class HandlerMixin:
    """Logging handler mixin class for acquiring and releasing a thread lock."""
    @staticmethod
    def lock(func: L) -> L:
        """A decorator for methods requiring acquire and release of a thread lock."""
        @wraps(func)
        def wrapper(self, *args, **kwargs) -> Any:
            self.acquire()
            try:
                return func(self, *args, **kwargs)
            finally:
                self.release()

        return cast(L, wrapper)

    @contextmanager
    def _locked(self) -> Generator[Any, Any, Any]:
        """A context manager for acquiring and releasing a thread lock."""
        self.acquire()
        try:
            yield
        finally:
            self.release()


class MsgpackHmacSocketHandler(logging.handlers.SocketHandler, HandlerMixin):
    """
    Subclass of `logging.handlers.SocketHandler` using MessagePack and hash-based message
    authentication codes (HMAC).
    """
    _supported_types = (str, int, float, bool, type(None), bytes, bytearray)

    def __init__(self, host: str, port: int) -> None:
        super().__init__(host, port)
        self.masked_keys: Dict[str, bytearray] = {}
        self.pack: types.BuiltinFunctionType = msgpack.Packer(use_bin_type=True).pack
        self._default_bytes: bytes = b""

    @HandlerMixin.lock
    def set_masked_key(self, task_id: str, masked_key: bytes) -> None:
        """Set a task ID and HMAC key into the cache."""
        self.masked_keys[task_id] = bytearray(masked_key)

    def pop_masked_key(self, task_id: str) -> None:
        """Pop a task ID and HMAC key from the cache."""
        with self._locked():
            masked_key = self.masked_keys.pop(task_id, None)
        self.zeroize(masked_key)

    @HandlerMixin.lock
    def clear_masked_keys(self) -> None:
        """Clear all task IDs and HMAC keys from the cache."""
        for task_id in list(self.masked_keys.keys()):
            masked_key = self.masked_keys.pop(task_id, None)
            self.zeroize(masked_key)
        self.masked_keys.clear()

    def zeroize(self, buffer: Optional[bytearray]) -> None:
        """Zeroize a bytearray in memory."""
        if isinstance(buffer, bytearray):
            buffer[:] = b"\x00" * len(buffer)

    def close(self) -> None:
        """Close the handler."""
        self.clear_masked_keys()
        super().close()

    def sanitize_record_arg(self, arg: Any) -> Any:
        """Sanitize a single element of log record `args` for MessagePack."""
        try:
            return arg if isinstance(arg, MsgpackHmacSocketHandler._supported_types) else str(arg)
        except Exception:
            return repr(arg)

    def sanitize_record_args(self, args: Any) -> Any:
        """Sanitize log record `args` for MessagePack."""
        if isinstance(args, dict):
            return args
        elif isinstance(args, list):
            return list(map(self.sanitize_record_arg, args))
        elif isinstance(args, tuple):
            return tuple(map(self.sanitize_record_arg, args))
        else:
            return self.sanitize_record_arg(args)

    def makePickle(self, record: logging.LogRecord) -> bytes:
        """Compress a logging record with MessagePack and a hash-based message authentication code (HMAC)."""
        record_dict = dict(
            msg=record.msg,
            args=self.sanitize_record_args(record.args),
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
            taskName=getattr(record, "taskName", None), # Added in Python-3.12
            protocol_name=record.protocol_name, # Set by `SetProtocolNameFilter` or `logging.LoggerAdapter`
            socket_address=record.socket_address, # Set by `SetSocketAddressFilter` or `logging.LoggerAdapter`
        )
        if record.exc_info:
            record_dict["exc_text"] = logging.Formatter().formatException(record.exc_info)
        if record.stack_info:
            record_dict["stack_info"] = record.stack_info

        packed_record = self.pack(record_dict)
        task_id = record.task_id # Set by `SetTaskIdFilter` or `logging.LoggerAdapter`
        with self._locked():
            masked_key = bytes(self.masked_keys.get(task_id, self._default_bytes))
        if masked_key == self._default_bytes:
            raise ValueError("`MsgpackHmacSocketHandler` could not get key from task ID.")
        frame = dict(task_id=task_id, packed_record=packed_record, version=1.0)
        packed_frame = self.pack(frame)
        signature = hmac_digest(masked_key, packed_frame)
        packet = dict(signature=signature, packed_frame=packed_frame)
        packed_packet = self.pack(packet)
        header = struct.pack(">L", len(packed_packet))

        return header + packed_packet


class MultiSocketHandler(logging.Handler, HandlerMixin):
    """
    Cache mutable dask worker logger handlers up to a maximum size, pruning least recently used (LRU)
    dask worker loggers first.
    """
    def __init__(self, logging_level=logging.NOTSET, maxsize=128) -> None:
        super().__init__()
        self.cache: OrderedDict[Tuple[str, int], MsgpackHmacSocketHandler] = collections.OrderedDict()
        self.maxsize: int = maxsize
        self.stdout_handler: logging.Handler = get_stdout_handler(logging_level=logging_level)

    def set_masked_key(self, socket_listener_address: Tuple[str, int], task_id: str, masked_key: bytes) -> None:
        """Set a masked key to handler cache."""
        host, port = socket_listener_address
        with self._locked():
            key, handler = self.get(host, port)
        handler.set_masked_key(task_id, masked_key)

    def pop_masked_key(self, socket_listener_address: Tuple[str, int], task_id: str) -> None:
        """Pop a masked key a handler cache."""
        host, port = socket_listener_address
        with self._locked():
            key, handler = self.get(host, port)
        handler.pop_masked_key(task_id)

    def setup_handler(self, host: str, port: int) -> MsgpackHmacSocketHandler:
        """Setup a `MsgpackHmacSocketHandler` instance."""
        handler = MsgpackHmacSocketHandler(host, port)
        handler.closeOnError = True

        return handler

    def get(self, host: str, port: int) -> Tuple[Tuple[str, int], MsgpackHmacSocketHandler]:
        """Set a key as most recently used, and return the key and value from the cache."""
        key = (host, port)
        with self._locked():
            handler = self.cache.pop(key, None) or self.setup_handler(host, port)
            self.cache[key] = handler # Most recently used
        self.maybe_prune()

        return key, handler

    def emit(self, record: logging.LogRecord) -> None:
        """Logging handler custom emit method override."""
        protocol_name = getattr(record, "protocol_name", None)
        socket_address = getattr(record, "socket_address", None)
        task_id = getattr(record, "task_id", None)
        if any(x in (None, "-", "-:0") for x in (protocol_name, socket_address, task_id)):
            # Log record was emitted from root logger on dask worker
            self.stdout_handler.handle(record)
            return
        host, port = split_socket_address(socket_address)
        with self._locked():
            key, handler = self.get(host, port)
        try:
            handler.handle(record)
        except Exception:
            # Drop broken socket so next log recreates it
            with self._locked():
                handler = self.cache.pop(key, None)
                if handler:
                    with suppress(Exception):
                        handler.close()

    def maybe_prune(self) -> None:
        """Prune the least recently used (LRU) items within the maximum size of the cache."""
        while len(self.cache) > self.maxsize:
            with self._locked():
                _, handler = self.cache.popitem(last=False)
            handler.flush()
            with suppress(Exception):
                handler.close()

    def purge_address(self, key: Tuple[str, int]) -> None:
        """Close and remove an item from the cache."""
        with self._locked():
            handler = self.cache.pop(key, None)
        if handler:
            handler.flush()
            with suppress(Exception):
                handler.close()

    @HandlerMixin.lock
    def purge_all(self) -> None:
        """Close and remove all items from the cache."""
        for handler in self.cache.values():
            handler.flush()
            with suppress(Exception):
                handler.close()
        self.cache.clear()

    def close(self) -> None:
        """Logging handler custom close method override."""
        self.purge_all()
        super().close()


def get_stdout_handler(logging_level=logging.NOTSET) -> logging.Handler:
    """Get a logging stream handler to `sys.stdout` for root logger records from dask workers."""
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging_level)
    formatter = logging.Formatter(
        ":".join(
            [
                "PyRosettaCluster_dask_worker_root",
                "%(asctime)s",
                "%(levelname)s",
                "%(name)s",
                " %(message)s",
            ]
        )
    )
    handler.setFormatter(formatter)

    return handler
