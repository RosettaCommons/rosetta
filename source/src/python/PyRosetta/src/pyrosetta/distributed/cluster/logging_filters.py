# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import logging

from typing import Tuple


class DefaultProtocolNameFilter(logging.Filter):
    """Maybe set protocol name of third-party root log records for logging socket listener formatter."""
    def __init__(self) -> None:
        super().__init__()

    def filter(self, record: logging.LogRecord) -> bool:
        if not hasattr(record, "protocol_name"):
            record.protocol_name = "-"
        return True


class SetProtocolNameFilter(logging.Filter):
    """Set protocol name for logging socket listener formatter."""
    def __init__(self, protocol_name: str) -> None:
        super().__init__()
        self.protocol_name = protocol_name

    def filter(self, record: logging.LogRecord) -> bool:
        record.protocol_name = self.protocol_name
        return True


class DefaultSocketAddressFilter(logging.Filter):
    """Maybe set socket address of third-party root log records for logging socket listener formatter."""
    def __init__(self) -> None:
        super().__init__()

    def filter(self, record: logging.LogRecord) -> bool:
        if not hasattr(record, "socket_address"):
            record.socket_address = format_socket_address(("-", 0))
        return True


class SetSocketAddressFilter(logging.Filter):
    """Set socket address for logging socket listener filter."""
    def __init__(self, socket_listener_address: Tuple[str, int]) -> None:
        super().__init__()
        self.socket_address = format_socket_address(socket_listener_address)

    def filter(self, record: logging.LogRecord) -> bool:
        record.socket_address = self.socket_address
        return True


class SocketAddressFilter(logging.Filter):
    """Filter log records for the logging socket listener address."""
    def __init__(self, socket_listener_address: Tuple[str, int]) -> None:
        super().__init__()
        self.socket_address = format_socket_address(socket_listener_address)

    def filter(self, record: logging.LogRecord) -> bool:
        return record.socket_address == self.socket_address


def format_socket_address(socket_listener_address: Tuple[str, int]) -> str:
    """Format a socket listener address for socket listener handler filters."""
    return ":".join(map(str, socket_listener_address))


def split_socket_address(socket_address: str) -> Tuple[str, int]:
    """Split a socket listener address for socket listener handlers."""
    host, port = tuple(s.strip() for s in socket_address.split(":"))

    return host, int(port)
