# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

import logging
import os
import warnings

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
        fh = logging.FileHandler(os.path.join(self.logs_path, "PyRosettaCluster.log",))
        fh.setFormatter(
            logging.Formatter(
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
        )
        logger.addHandler(fh)

    def _close_logger(self) -> None:
        """Close the logger for the client instance."""

        logger = logging.getLogger()
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
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
        logging_file,
        logging_level,
        DATETIME_FORMAT,
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

        logs_path = os.path.dirname(logging_file)
        if not os.path.isdir(logs_path):
            warnings.warn(
                f"Creating logs directory in 'LoggingSupport.setup_target_logging': {logs_path}",
                RuntimeWarning,
                stacklevel=2,
            )
            os.makedirs(logs_path, exist_ok=True)

        pid = os.getpid()
        logging_file_split = list(os.path.splitext(logging_file))
        logging_file_split.insert(-1, f"-{pid}")
        logging_file = "".join(logging_file_split)
        print("Setting up logging file on dask worker:", logging_file)
        fh = logging.FileHandler(logging_file, delay=True)
        fh.setFormatter(
            logging.Formatter(
                ":".join(
                    [
                        "%(levelname)s",
                        str(protocol.__name__),
                        "%(name)s",
                        "%(asctime)s",
                        " %(message)s",
                    ]
                )
            )
        )
        logger.addHandler(fh)

        result = func(
            protocol,
            compressed_packed_pose,
            compressed_kwargs,
            q,
            logging_file,
            logging_level,
            DATETIME_FORMAT,
            ignore_errors,
            protocols_key,
            decoy_ids,
            compression,
            client_residue_type_set,
            client_repr,
            **pyrosetta_init_kwargs,
        )

        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            with suppress(Exception):
                handler.close()

        return result

    return cast(L, wrapper)
