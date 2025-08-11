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

from contextlib import suppress
from typing import (
    Any,
    Dict,
    Optional,
    OrderedDict,
    Tuple,
)

from pyrosetta.distributed.cluster.hkdf import hmac_digest
from pyrosetta.distributed.cluster.logging_filters import split_socket_address


class MsgpackHmacSocketHandler(logging.handlers.SocketHandler):
    """
    Subclass of `logging.handlers.SocketHandler` using MessagePack and hash-based message
    authentication codes (HMAC).
    """
    _supported_types = (str, int, float, bool, type(None), bytes, bytearray)

    def __init__(self, host: str, port: int) -> None:
        super().__init__(host, port)
        self.masked_keys: Dict[str, bytearray] = {}

    def set_masked_key(self, task_id: str, masked_key: bytes) -> None:
        """Set a task ID and HMAC key into the cache."""
        self.acquire()
        try:
            self.masked_keys[task_id] = bytearray(masked_key)
        finally:
            self.release()

    def pop_masked_key(self, task_id: str) -> None:
        """Pop a task ID and HMAC key from the cache."""
        self.acquire()
        try:
            masked_key = self.masked_keys.pop(task_id, None)
            self.zeroize(masked_key)
        finally:
            self.release()

    def clear_masked_keys(self) -> None:
        """Clear all task IDs and HMAC keys from the cache."""
        self.acquire()
        try:
            for task_id in list(self.masked_keys.keys()):
                masked_key = self.masked_keys.pop(task_id, None)
                self.zeroize(masked_key)
            self.masked_keys.clear()
        finally:
            self.release()

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

        packed_record = msgpack.packb(record_dict, use_bin_type=True)
        task_id = record.task_id # Set by `SetTaskIdFilter` or `logging.LoggerAdapter`
        self.acquire()
        try: # Acquire while getting masked key and HMAC digest
            masked_key = self.masked_keys.get(task_id, None)
            if masked_key is None:
                raise ValueError("MsgpackHmacSocketHandler could not get key from task ID.")
            frame = dict(task_id=task_id, packed_record=packed_record, version=1.0)
            packed_frame = msgpack.packb(frame, use_bin_type=True)
            signature = hmac_digest(masked_key, packed_frame)
        except ValueError as ex:
            raise ValueError(ex)
        finally:
            self.release()
        packet = dict(signature=signature, packed_frame=packed_frame)
        packed_packet = msgpack.packb(packet, use_bin_type=True)
        header = struct.pack(">L", len(packed_packet))

        return header + packed_packet


class MultiSocketHandler(logging.Handler):
    """
    Cache mutable dask worker logger handlers up to a maximum size, pruning least recently used (LRU)
    dask worker loggers first.
    """
    def __init__(self, maxsize=128) -> None:
        super().__init__()
        self.cache: OrderedDict[Tuple[str, int], MsgpackHmacSocketHandler] = collections.OrderedDict()
        self.maxsize: int = maxsize
        self.stdout_handler: logging.Handler = get_stdout_handler()

    def set_masked_key(self, socket_listener_address: Tuple[str, int], task_id: str, masked_key: bytes) -> None:
        """Set a masked key to handler cache."""
        host, port = socket_listener_address
        self.acquire()
        try:
            key, handler = self.get(host, port)
            handler.set_masked_key(task_id, masked_key)
        finally:
            self.release()

    def pop_masked_key(self, socket_listener_address: Tuple[str, int], task_id: str) -> None:
        """Pop a masked key a handler cache."""
        host, port = socket_listener_address
        self.acquire()
        try:
            key, handler = self.get(host, port)
            handler.pop_masked_key(task_id)
        finally:
            self.release()

    def setup_handler(self, host: str, port: int) -> MsgpackHmacSocketHandler:
        """Setup a `MsgpackHmacSocketHandler` instance."""
        handler = MsgpackHmacSocketHandler(host, port)
        handler.closeOnError = True

        return handler

    def get(self, host: str, port: int) -> Tuple[Tuple[str, int], MsgpackHmacSocketHandler]:
        """Set a key as most recently used, and return the key and value from the cache."""
        key = (host, port)
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
        self.acquire()
        try:
            key, handler = self.get(host, port)
            try:
                handler.emit(record)
            except Exception:
                # Drop broken socket so next log recreates it
                with suppress(Exception):
                    handler.close()
                self.cache.pop(key, None)
        finally:
            self.release()

    def maybe_prune(self) -> None:
        """Prune the least recently used (LRU) items within the maximum size of the cache."""
        self.acquire()
        try:
            while len(self.cache) > self.maxsize:
                _, handler = self.cache.popitem(last=False)
                handler.flush()
                with suppress(Exception):
                    handler.close()
        finally:
            self.release()

    def purge_address(self, key: Tuple[str, int]) -> None:
        """Close and remove an item from the cache."""
        self.acquire()
        try:
            handler = self.cache.pop(key, None)
            if handler:
                handler.flush()
                with suppress(Exception):
                    handler.close()
        finally:
            self.release()

    def purge_all(self) -> None:
        """Close and remove all items from the cache."""
        self.acquire()
        try:
            for handler in self.cache.values():
                handler.flush()
                with suppress(Exception):
                    handler.close()
            self.cache.clear()
        finally:
            self.release()

    def close(self) -> None:
        """Logging handler custom close method override."""
        self.purge_all()
        super().close()


def get_stdout_handler() -> logging.Handler:
    """Get a logging stream handler to `sys.stdout` for root logger records from dask workers."""
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.NOTSET)
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
