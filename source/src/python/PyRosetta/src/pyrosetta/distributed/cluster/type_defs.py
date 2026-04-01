# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import sys

if sys.version_info[:2] < (3, 9):
    from typing import Callable, Sequence
else:
    from collections.abc import Callable, Sequence

from pyrosetta import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    Dict,
    Generator,
    List,
    Optional,
    Tuple,
    TypeVar,
    Union,
)


class Pose(Pose):
    """
    Subclass of `Pose` object to coax type checkers into treating `Pose` as a proper class type instead of type
    `Any` due to the underlying C++ bindings being opaque to type checkers.
    """
    def __init__(self) -> None:
        super().__init__()


# User-defined PyRosetta protocol type abstractions:

PyRosettaProtocolResult = Union[PackedPose, Pose, Dict[str, Any], None]
"""An individual output result from a user-defined PyRosetta protocol."""

PyRosettaProtocolGenerator = Generator[PyRosettaProtocolResult, None, None]
"""A generator yielding `PyRosettaProtocolResult` objects."""

PyRosettaProtocolResults = Union[
    PyRosettaProtocolResult,
    Sequence[PyRosettaProtocolResult],
    PyRosettaProtocolGenerator,
]
"""Collective output results from a user-defined PyRosetta protocol."""

# Generic type abstractions:

T = TypeVar("T")
"""An arbitrary type."""

ListOrTuple = Union[List[T], Tuple[T, ...]]
"""
A container that must be either a `list` or `tuple` object, where all elements follow the same declared type
pattern `T`.
"""

FloatOrInt = Union[float, int]
"""A `float` or `int` object."""

PoseOrPackedPose = Union[Pose, PackedPose]
"""A `Pose` or `PackedPose` object."""

# Internal method type abstractions:

PyRosettaProtocol = Callable[..., PyRosettaProtocolResults]
"""A callable user-defined PyRosetta protocol."""

PyRosettaProtocols = Sequence[PyRosettaProtocol]
"""A sequence of callable user-defined PyRosetta protocols."""

TaskResource = Optional[Dict[str, FloatOrInt]]
"""An optional Dask resource constraint for an individual task."""

TaskChainResource = Optional[ListOrTuple[TaskResource]]
"""An optional container of Dask resource constraints for a task chain."""

TaskChainClientIndices = Optional[ListOrTuple[int]]
"""An optional container of Dask client indices for a task chain."""

TaskChainPriorities = Optional[ListOrTuple[int]]
"""An optional container of priorities for a task chain."""

TaskChainRetries = Optional[Union[int, ListOrTuple[int]]]
"""An optional `int` object or container of retries for a task chain."""
