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
import hashlib
import hmac
import logging
import struct
import sys

from contextlib import suppress
from typing import (
    Optional,
    OrderedDict,
    Tuple,
)

from pyrosetta.distributed.cluster.logging_filters import split_socket_address


class MsgpackHmacSocketHandler(logging.handlers.SocketHandler):
    """
    Subclass of `logging.handlers.SocketHandler` using MessagePack and hash-based message
    authentication codes (HMAC).
    """
    def __init__(self, host: str, port: int, masked_key: Optional[bytes] = None) -> None:
        super().__init__(host, port)
        self.masked_key = masked_key

    def set_masked_key(self, masked_key: bytes) -> None:
        self.masked_key = masked_key

    def makePickle(self, record: logging.LogRecord) -> bytes:
        """Compress a logging record with MessagePack and a hash-based message authentication code (HMAC)."""
        if not isinstance(self.masked_key, bytes):
            raise TypeError(f"Cannot sign HMAC with key of type: {type(self.masked_key)}")

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
            taskName=record.taskName if sys.version_info[:2] >= (3, 12) else None,
            protocol_name=record.protocol_name, # Set by `SetProtocolNameFilter` or `logging.LoggerAdapter`
            socket_address=record.socket_address, # Set by `SetSocketAddressFilter` or `logging.LoggerAdapter`
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


class MultiSocketHandler(logging.Handler):
    """
    Cache mutable dask worker logger handlers up to a maximum size, pruning least recently used (LRU)
    dask worker loggers first.
    """
    def __init__(self, maxsize=128):
        super().__init__()
        self.cache: OrderedDict[Tuple[str, int], MsgpackHmacSocketHandler] = collections.OrderedDict()
        self.maxsize: int = maxsize

    def set_masked_key(self, socket_listener_address, masked_key) -> None:
        host, port = socket_listener_address
        self.acquire()
        try:
            key, handler = self.get(host, port)
            handler.set_masked_key(masked_key)
        finally:
            self.release()

    def setup_handler(self, host: str, port: int) -> MsgpackHmacSocketHandler:
        handler = MsgpackHmacSocketHandler(host, port)
        handler.closeOnError = True

        return handler

    def get(self, host: str, port: int) -> Tuple[Tuple[str, int], MsgpackHmacSocketHandler]:
        key = (host, port)
        handler = self.cache.pop(key, None) or self.setup_handler(host, port)
        self.cache[key] = handler # Most recently used
        self.maybe_prune()

        return key, handler

    def emit(self, record: logging.LogRecord) -> None:
        """Logging handler custom emit method override."""
        socket_address = getattr(record, "socket_address", None)
        if not socket_address:
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
