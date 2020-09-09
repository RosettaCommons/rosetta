# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"
__email__ = "klima.jason@gmail.com"

import logging
import os

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
        """Open the logger for the master instance."""

        logger = logging.getLogger()
        logger.setLevel(self.logging_level)
        logger.handlers = []
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
        """Close the logger for the master instance."""

        logger = logging.getLogger()
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)


def setup_target_logging(func: L) -> L:
    """Support logging within the spawned thread."""

    @wraps(func)
    def wrapper(
        protocol,
        pose,
        q,
        logging_file,
        logging_level,
        DATETIME_FORMAT,
        ignore_errors,
        master_residue_type_set,
        **kwargs,
    ):
        """Wrapper function to setup_target_logging."""

        logging.basicConfig(
            filename=logging_file,
            level=logging_level,
            format=":".join(
                [
                    "%(levelname)s",
                    str(protocol.__name__),
                    "%(name)s",
                    "%(asctime)s",
                    " %(message)s",
                ]
            ),
        )

        return func(
            protocol,
            pose,
            q,
            logging_file,
            logging_level,
            DATETIME_FORMAT,
            ignore_errors,
            master_residue_type_set,
            **kwargs,
        )

    return cast(L, wrapper)
