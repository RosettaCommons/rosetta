# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.base' requires the "
        + "third-party package 'toolz' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import pyrosetta
import pyrosetta.distributed

from datetime import datetime
from functools import wraps
from pyrosetta.distributed.cluster.converters import _parse_protocols
from pyrosetta.distributed.cluster.initialization import (
    _get_residue_type_set_name3 as _get_residue_type_set,
)
from pyrosetta.distributed.cluster.validators import (
    _validate_clients_indices,
    _validate_protocols_seeds_decoy_ids,
    _validate_resources,
)
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
        task: Dict[Any, Any],
    ) -> Tuple[bytes, Dict[str, Any]]:
        """Setup the kwargs for the initial task."""

        kwargs = {
            self.protocols_key: protocols.copy(),
            "PyRosettaCluster_output_path": self.output_path,
            "PyRosettaCluster_logging_file": self.logging_file,
            "PyRosettaCluster_task": self.serializer.deepcopy_kwargs(task),
            **self._setup_seed(self.serializer.deepcopy_kwargs(task), seed),
        }
        pyrosetta_init_kwargs = self._setup_pyrosetta_init_kwargs(kwargs)
        compressed_kwargs = self.serializer.compress_kwargs(kwargs)

        return compressed_kwargs, pyrosetta_init_kwargs

    def _setup_pyrosetta_init_kwargs(self, kwargs: Dict[Any, Any]) -> Dict[str, Any]:
        pyrosetta_init_kwargs = toolz.dicttoolz.keyfilter(
            lambda k: k in self.pyrosetta_init_args,
            kwargs,
        )

        return pyrosetta_init_kwargs
    
    def _get_clients_index(
            self, clients_indices: List[int], protocols: List[Callable[..., Any]]
        ) -> int:
        """Return the clients index for the current protocol."""
        if clients_indices is None:
            return 0
        else:
            _protocols_index = len(clients_indices) - len(protocols)
            return clients_indices[_protocols_index]
    
    def _get_resource(
            self, resources: List[Dict[Any, Any]], protocols: List[Callable[..., Any]]
        ) -> Optional[Dict[Any, Any]]:
        """Return the resource for the current protocol."""
        if resources is None:
            return None
        else:
            _protocols_index = len(resources) - len(protocols)
            return resources[_protocols_index]

    def _setup_kwargs(
        self, kwargs: Dict[Any, Any], clients_indices: List[int], resources: Optional[Dict[Any, Any]],
    ) -> Tuple[bytes, Dict[str, Any], Callable[..., Any], int, Optional[Dict[Any, Any]]]:
        """Setup the kwargs for the subsequent tasks."""
        clients_index = self._get_clients_index(clients_indices, kwargs[self.protocols_key])
        resource = self._get_resource(resources, kwargs[self.protocols_key])
        _protocols, protocol, seed = self._get_task_state(kwargs[self.protocols_key])
        kwargs[self.protocols_key] = _protocols
        kwargs = self._setup_seed(kwargs, seed)
        pyrosetta_init_kwargs = self._setup_pyrosetta_init_kwargs(kwargs)
        compressed_kwargs = self.serializer.compress_kwargs(kwargs)

        return compressed_kwargs, pyrosetta_init_kwargs, protocol, clients_index, resource

    def _setup_seed(self, kwargs: Dict[Any, Any], seed: Optional[str]) -> Dict[Any, Any]:
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
        self, args: Tuple[Any, ...], protocols: Any, clients_indices: Any, resources: Any,
    ) -> Tuple[List[Callable[..., Any]], Callable[..., Any], Optional[str], int, Optional[Dict[Any, Any]]]:
        """Parse, validate, and setup the user-provided PyRosetta protocol(s)."""

        _protocols = _parse_protocols(args) + _parse_protocols(protocols)
        _clients_dict_keys = list(self.clients_dict.keys())
        _validate_clients_indices(clients_indices, _protocols, _clients_dict_keys)
        _validate_resources(resources, _protocols)
        _clients_index = self._get_clients_index(clients_indices, _protocols)
        _resource = self._get_resource(resources, _protocols)

        _protocols, _protocol, _seed = self._get_task_state(
            _validate_protocols_seeds_decoy_ids(_protocols, self.seeds, self.decoy_ids)
        )

        return _protocols, _protocol, _seed, _clients_index, _resource


def capture_task_metadata(func: M) -> M:
    """Capture a task's metadata as kwargs."""

    @wraps(func)
    def wrapper(
        protocol,
        compressed_packed_pose,
        DATETIME_FORMAT,
        ignore_errors,
        protocols_key,
        decoy_ids,
        serializer,
        **kwargs,
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

        return func(
            protocol,
            compressed_packed_pose,
            DATETIME_FORMAT,
            ignore_errors,
            protocols_key,
            decoy_ids,
            serializer,
            **kwargs,
        )

    return cast(M, wrapper)
