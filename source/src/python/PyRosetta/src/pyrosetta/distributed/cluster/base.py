# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"
__email__ = "klima.jason@gmail.com"

import copy
import logging
import pyrosetta
import pyrosetta.distributed

from datetime import datetime
from functools import wraps
from pyrosetta.distributed.cluster.converters import _parse_protocols
from pyrosetta.distributed.cluster.initialization import (
    _get_residue_type_set_name3 as _get_residue_type_set,
)
from pyrosetta.distributed.cluster.validators import _validate_protocols_seeds_decoy_ids
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    Callable,
    Dict,
    Generic,
    List,
    Optional,
    Sized,
    Tuple,
    TypeVar,
    Union,
    cast,
)


G = TypeVar("G")
M = TypeVar("M", bound=Callable[..., Any])


class TaskBase(Generic[G]):
    """Task objects underpinning PyRosettaCluster."""

    def _get_seed(self, protocols: Sized) -> Optional[str]:
        """Get the seed for the input user-provided PyRosetta protocol."""

        if self.seeds:
            seed_index = (len(self.seeds) - len(protocols)) - 1
            seed = self.seeds[seed_index]
        else:
            seed = None

        return seed

    def _get_task_state(
        self, protocols: List[Callable[..., Any]]
    ) -> Tuple[List[Callable[..., Any]], Callable[..., Any], Optional[str]]:
        """
        Given the current state of protocols, returns a tuple of the updated
        state of protocols and current protocol and seed.
        """

        protocol = protocols.pop(0)
        seed = self._get_seed(protocols)

        return protocols, protocol, seed

    def _setup_initial_kwargs(
        self,
        protocols: List[Callable[..., Any]],
        seed: Optional[str],
        kwargs: Dict[Any, Any],
    ) -> Dict[str, Union[Dict[Any, Any], List[Callable[..., Any]], str]]:
        """Setup the kwargs for the initial task."""

        return {
            self.protocols_key: protocols.copy(),
            "PyRosettaCluster_logging_file": self.logging_file,
            "PyRosettaCluster_task": copy.deepcopy(kwargs),
            **self._setup_seed(copy.deepcopy(kwargs), seed),
        }

    def _setup_kwargs(
        self, kwargs: Dict[Any, Any]
    ) -> Tuple[Dict[Any, Any], Callable[..., Any]]:
        """Setup the kwargs for the subsequent tasks."""

        _protocols, protocol, seed = self._get_task_state(kwargs[self.protocols_key])
        kwargs[self.protocols_key] = _protocols
        kwargs = self._setup_seed(kwargs, seed)

        return kwargs, protocol

    def _setup_seed(
        self, kwargs: Dict[Any, Any], seed: Optional[str]
    ) -> Dict[Any, Any]:
        """
        Setup the 'options' or 'extra_options' task kwargs with the `-run:jran`
        PyRosetta command line flag.
        """

        if seed:
            jran = f" -run:constant_seed 1 -run:jran {seed}"
            split_str = " -run:constant_seed 1 -run:jran "
            if "extra_options" in kwargs:
                flags = pyrosetta.distributed._normflags(kwargs["extra_options"])
                if split_str in flags:
                    flags = flags.split(split_str)[0]
                kwargs["extra_options"] = flags + jran
            elif "options" in kwargs:
                flags = pyrosetta.distributed._normflags(kwargs["options"])
                if split_str in flags:
                    flags = flags.split(split_str)[0]
                kwargs["options"] = (
                    pyrosetta.distributed._normflags(kwargs["options"]) + jran
                )
            else:
                kwargs["extra_options"] = "-out:levels all:warning" + jran

        return kwargs

    def _setup_protocols_protocol_seed(
        self, args: Tuple[Any, ...], protocols: Any
    ) -> Tuple[List[Callable[..., Any]], Callable[..., Any], Optional[str]]:
        """Parse, validate, and setup the user-provided PyRosetta protocol(s)."""

        _protocols = _parse_protocols(args) + _parse_protocols(protocols)

        return self._get_task_state(
            _validate_protocols_seeds_decoy_ids(_protocols, self.seeds, self.decoy_ids)
        )


def _get_decoy_id(protocols: Sized, decoy_ids: List[int]) -> Optional[int]:
    """Get the decoy number given the user-provided PyRosetta protocols."""

    if decoy_ids:
        decoy_id_index = (len(decoy_ids) - len(protocols)) - 1
        decoy_id = decoy_ids[decoy_id_index]
    else:
        decoy_id = None

    return decoy_id


def _get_packed_pose_kwargs_pairs_list(
    results: List[PackedPose],
    kwargs: Dict[Any, Any],
    protocol_name: str,
    protocols_key: str,
    decoy_ids: List[int],
) -> List[Tuple[PackedPose, Dict[Any, Any]]]:
    """Parse PackedPose and kwargs objects into a list of tuples."""

    decoy_id = _get_decoy_id(kwargs[protocols_key], decoy_ids)
    packed_pose_kwargs_pairs_list = []
    for i, packed_pose in enumerate(results):
        if (decoy_id != None) and (i != decoy_id):
            logging.info(
                "Discarding a returned decoy because it does "
                + "not match the user-provided decoy_ids."
            )
            continue
        task_kwargs = copy.deepcopy(kwargs)
        if "PyRosettaCluster_decoy_ids" not in task_kwargs:
            task_kwargs["PyRosettaCluster_decoy_ids"] = []
        task_kwargs["PyRosettaCluster_decoy_ids"].append((protocol_name, i))
        packed_pose_kwargs_pairs_list.append((packed_pose, task_kwargs))
    if decoy_id:
        assert (
            len(packed_pose_kwargs_pairs_list) == 1
        ), "When specifying decoy_ids, there may only be one decoy_id per protocol."

    return packed_pose_kwargs_pairs_list


def capture_task_metadata(func: M) -> M:
    """Capture a task's metadata as kwargs."""

    @wraps(func)
    def wrapper(
        protocol, pose, DATETIME_FORMAT, ignore_errors, **kwargs,
    ):
        """Wrapper function to capture_task_metadata."""

        protocol_name = str(protocol.__name__)
        kwargs["PyRosettaCluster_protocol_name"] = protocol_name
        if not "PyRosettaCluster_protocols" in kwargs:
            kwargs["PyRosettaCluster_protocols"] = []
        kwargs["PyRosettaCluster_protocols"].append(protocol_name)
        kwargs["PyRosettaCluster_protocol_number"] = (
            len(kwargs["PyRosettaCluster_protocols"]) - 1
        )  # Pythonic indices
        if not "PyRosettaCluster_datetime_start" in kwargs:
            kwargs["PyRosettaCluster_datetime_start"] = datetime.now().strftime(
                DATETIME_FORMAT
            )
        if not "PyRosettaCluster_seeds" in kwargs:
            kwargs["PyRosettaCluster_seeds"] = []
        kwargs["PyRosettaCluster_seeds"].append(
            (protocol_name, str(pyrosetta.rosetta.numeric.random.rg().get_seed()))
        )

        return func(protocol, pose, DATETIME_FORMAT, ignore_errors, **kwargs,)

    return cast(M, wrapper)
