# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import sys

if sys.version_info[:2] < (3, 9):
    from typing import (
        AbstractSet,
        Callable,
        Deque,
    )
else:
    from collections import deque as Deque
    from collections.abc import (
        Callable,
        Set as AbstractSet,
    )

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

# Generic type abstractions:

T = TypeVar("T")
"""An arbitrary type variable."""

CallableType = TypeVar("CallableType", bound=Callable[..., Any])
"""An arbitrary callable type variable."""

ListOrTuple = Union[List[T], Tuple[T, ...]]
"""
A container that must be either a `list` or `tuple` object, where all elements follow the same declared type
pattern `T`.
"""

FloatOrInt = Union[float, int]
"""A `float` or `int` object."""

PoseOrPackedPose = Union[Pose, PackedPose]
"""A `Pose` or `PackedPose` object."""

# User-defined PyRosetta protocol type abstractions:

PyRosettaProtocolResult = Union[PoseOrPackedPose, Dict[str, Any], None]
"""An individual output result from a user-defined PyRosetta protocol."""

PyRosettaProtocolGenerator = Generator[PyRosettaProtocolResult, None, None]
"""A generator yielding `PyRosettaProtocolResult` objects."""

PyRosettaProtocolResults = Union[
    PyRosettaProtocolResult,
    ListOrTuple[PyRosettaProtocolResult],
    PyRosettaProtocolGenerator,
]
"""Collective output results from a user-defined PyRosetta protocol."""

PyRosettaProtocol = Callable[..., PyRosettaProtocolResults]
"""A callable user-defined PyRosetta protocol."""

# Internal method type abstractions:

PyRosettaProtocols = List[PyRosettaProtocol]
"""A `list` object of callable user-defined PyRosetta protocols."""

TaskResource = Optional[Dict[str, FloatOrInt]]
"""An optional Dask resource constraint for an individual task."""

TaskChainResources = Optional[ListOrTuple[TaskResource]]
"""An optional container of Dask resource constraints for a task chain."""

TaskChainClientIndices = Optional[ListOrTuple[int]]
"""An optional container of Dask client indices for a task chain."""

TaskChainPriorities = Optional[ListOrTuple[int]]
"""An optional container of priorities for a task chain."""

TaskChainRetries = Optional[Union[int, ListOrTuple[int]]]
"""An optional `int` object or container of retries for a task chain."""
