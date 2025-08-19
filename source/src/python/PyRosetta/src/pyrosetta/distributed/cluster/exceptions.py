# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    import billiard
    from billiard import WorkerLostError
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.exceptions' requires the "
        + "third-party package 'billiard' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/billiard/\n"
    )
    raise

import logging
import traceback

from functools import wraps
from pyrosetta.distributed.packed_pose.core import PackedPose
from queue import Empty
from typing import (
    Any,
    Callable,
    Dict,
    List,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
    cast,
)


Q = TypeVar("Q", bound=billiard.Queue)
P = TypeVar("P", bound=billiard.context.Process)
T = TypeVar("T", bound=Callable[..., Any])


class InputError(TypeError):
    """
    Exception raised for PyRosettaCluster attribute input errors for `str` and `int` types.
    """

    def __init__(self, obj: Any, attribute: str) -> NoReturn:
        super().__init__(
            f"PyRosettaCluster '{attribute}' attribute must be of type `int` or `str`, "
            + f"or an iterable containing `int` or `str` types. Received {obj} of type `{type(obj)}`."
        )


class InputFileError(TypeError):
    """
    Exception raised for PyRosettaCluster  errors for `str` and `int` types.
    """

    def __init__(self, obj: Any) -> NoReturn:
        super().__init__(
            "The `input_file` argument parameter must be of type `str`, not of type {0}.".format(
                type(obj)
            )
        )


class OutputError(TypeError):
    """Exception raised for worker output errors."""

    def __init__(self, obj: Any) -> NoReturn:
        super().__init__(
            " ".join(
                "Returned object(s) should be an instance of `NoneType`, `Pose`, `PackedPose`, or `dict`; or \
                an iterable containing `NoneType`, `Pose`, or `PackedPose` objects, and optionally one object of type `dict`. \
                The output object obtained was `{0}` of type `{1}`.".format(
                    obj, type(obj)
                ).split()
            )
        )


class WorkerError(WorkerLostError):
    """Exception raised for worker errors."""

    def __init__(self, protocol_name: str) -> NoReturn:
        super().__init__(WorkerError._msg(protocol_name))

    @staticmethod
    def _msg(protocol_name: str) -> str:
        return (
            "Worker thread killed due to an error or segmentation fault encountered "
            + f"in the user-provided PyRosetta protocol '{protocol_name}'."
        )

    @staticmethod
    def _ignore_errors_msg(protocol_name: str) -> str:
        return (
            WorkerError._msg(protocol_name)
            + "Ignoring error because `ignore_errors` is enabled!"
        )

    @staticmethod
    def _get_msg(protocol_name: str, ignore_errors: bool) -> str:
        if ignore_errors:
            return WorkerError._ignore_errors_msg(protocol_name)
        else:
            return WorkerError._msg(protocol_name)


def trace_protocol_exceptions(func: T) -> Union[T, NoReturn]:
    """Trace exceptions in user-provided PyRosetta protocols."""
    @wraps(func)
    def wrapper(
        packed_pose: PackedPose,
        protocol: Callable[..., Any],
        ignore_errors: bool,
        **kwargs: Dict[Any, Any],
    ) -> Union[Any, NoReturn]:
        protocol_name = protocol.__name__
        try:
            result = func(packed_pose, protocol, ignore_errors, **kwargs)
        except:
            logging.error(
                traceback.format_exc()
                + WorkerError._get_msg(protocol_name, ignore_errors)
            )
            if ignore_errors:
                # Return a `NoneType` object to be converted to an empty `PackedPose` object
                # when a non-system-exiting Python exception is raised
                result = None
            else:
                raise

        return result

    return cast(T, wrapper)


def trace_subprocess_exceptions(func: T) -> Union[T, NoReturn]:
    """Trace exceptions in billiard subprocesses."""
    @wraps(func)
    def wrapper(
        q: Q,
        p: P,
        compressed_kwargs: bytes,
        protocol_name: str,
        timeout: Union[float, int],
        ignore_errors: bool,
    ) -> Union[List[Tuple[Optional[bytes], bytes]], NoReturn]:
        while True:
            try:
                _results = func(
                    q, p, compressed_kwargs, protocol_name, timeout, ignore_errors
                )
                break
            except Empty:
                if not p.is_alive():
                    if ignore_errors:
                        logging.error(WorkerError._ignore_errors_msg(protocol_name))
                        _results = func(
                            q, p, compressed_kwargs, protocol_name, timeout, ignore_errors
                        )
                        break
                    else:
                        raise WorkerError(protocol_name)

        return _results

    return cast(T, wrapper)
