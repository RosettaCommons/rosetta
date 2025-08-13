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
        "Importing 'pyrosetta.distributed.cluster.logging_listeners' requires the "
        + "third-party package 'msgpack' as a dependency!\n"
        + "Please install the package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/msgpack/\n"
    )
    raise

import logging
import os
import select
import socketserver
import struct
import threading
import traceback
import warnings

from functools import partial
from typing import (
    Any,
    Dict,
    Union,
)

from pyrosetta.distributed.cluster.hkdf import compare_digest, derive_task_key, hmac_digest
from pyrosetta.distributed.cluster.logging_filters import (
    SocketAddressFilter,
    split_socket_address,
)


class LogRecordRequestHandler(socketserver.StreamRequestHandler):
    """
    Handler for a streaming logging request modified from logging cookbook recipe:
    https://docs.python.org/3/howto/logging-cookbook.html#sending-and-receiving-logging-events-across-a-network
    """
    def setup(self) -> None:
        """Setup socket server."""
        super().setup()
        self.unpack: partial = partial(msgpack.unpackb, raw=False)

    def handle(self) -> None:
        """Handle raw logging message and make log record."""
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
        """
        Log record decompress method override using MessagePack and hash-based message authentication codes (HMAC).
        """
        packet = self.unpack(msg)
        signature = packet["signature"]
        packed_frame = packet["packed_frame"]
        frame = self.unpack(packed_frame)
        task_id = frame["task_id"]
        packed_record = frame["packed_record"]
        masked_key = derive_task_key(self.server.passkey, task_id)
        required_signature = hmac_digest(masked_key, packed_frame)
        if not compare_digest(required_signature, signature):
            raise ValueError("Logging socket listener received a bad hash-based message authentication code!")
        record = self.unpack(packed_record)
        # Update record.args for positional % formatting
        args = record.get("args", ())
        if isinstance(args, list):
            record["args"] = tuple(args)

        return record


class MaskedBytes(bytes):
    """A `bytes` subclass to mask its contents if the `repr` method is called."""
    def __new__(cls, value: bytes) -> bytes:
        return super().__new__(cls, value)

    def __repr__(self) -> str:
        return "#"


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
        _host, _port = split_socket_address(logging_address)
        super().__init__((_host, _port), LogRecordRequestHandler)
        self.socket_listener_address = self.socket.getsockname()
        self.handler = self.setup_handler(logging_file, logging_level)
        self.passkey = MaskedBytes(os.urandom(32))
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
        """Start logging socket listener."""
        self._thread = threading.Thread(target=self.serve_forever, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        """Stop logging socket listener."""
        self.shutdown()
        self.server_close()
        if self._thread:
            self._thread.join()

    def serve_until_stopped(self) -> None:
        """Run logging socket listener until aborted."""
        abort = 0
        while not abort:
            _read_ready, _, _ = select.select([self.socket.fileno()], [], [], self.timeout)
            if _read_ready:
                self.handle_request()
            abort = self.abort
