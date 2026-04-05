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
        + "Please install this package into your virtual environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import logging
import pyrosetta
import pyrosetta.distributed

from datetime import datetime
from functools import wraps
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Sized,
    Tuple,
    TypeVar,
    cast,
)

from pyrosetta.distributed.cluster.config import __dask_version__
from pyrosetta.distributed.cluster.converters import (
    _parse_protocols,
    _version_tuple_to_str,
)
from pyrosetta.distributed.cluster.initialization import (
    _get_norm_task_options,
    _get_residue_type_set_name3 as _get_residue_type_set,
)
from pyrosetta.distributed.cluster.serialization import Serialization
from pyrosetta.distributed.cluster.type_defs import (
    CallableType,
    PyRosettaProtocol,
    PyRosettaProtocols,
    TaskChainClientIndices,
    TaskChainPriorities,
    TaskChainResources,
    TaskChainRetries,
    TaskResource,
)
from pyrosetta.distributed.cluster.validators import (
    _validate_clients_indices,
    _validate_priorities,
    _validate_protocols_seeds_decoy_ids,
    _validate_resources,
    _validate_retries,
)


class TaskBase:
    """Task objects underpinning `PyRosettaCluster`."""

    def _get_seed(self, protocols: Sized) -> Optional[str]:
        """Get the PyRosetta RNG seed for the input user-defined PyRosetta protocol."""

        if self.seeds:
            seed_index = (len(self.seeds) - len(protocols)) - 1
            seed = self.seeds[seed_index]
        else:
            seed = None

        return seed

    def _get_task_state(
        self, protocols: PyRosettaProtocols
    ) -> Tuple[PyRosettaProtocols, PyRosettaProtocol, Optional[str]]:
        """
        Given the current state of the PyRosetta protocols, return a `tuple` object of the updated state of the
        PyRosetta protocols, the current PyRosetta protocol, and curent PyRosetta RNG seed.
        """

        protocol = protocols.pop(0)
        seed = self._get_seed(protocols)

        return protocols, protocol, seed

    def _setup_initial_kwargs(
        self,
        protocols: PyRosettaProtocols,
        seed: Optional[str],
        task: Dict[str, Any],
    ) -> Tuple[bytes, Dict[str, Any]]:
        """Setup the keyword arguments for the initial user-defined task dictionary."""

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

    def _setup_pyrosetta_init_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """Setup the `pyrosetta.init` function keyword arguments."""

        pyrosetta_init_kwargs = toolz.dicttoolz.keyfilter(
            lambda k: k in self.pyrosetta_init_args,
            kwargs,
        )

        return pyrosetta_init_kwargs

    def _get_clients_index(
        self, clients_indices: TaskChainClientIndices, protocols: Sized
    ) -> int:
        """Return the clients index for the current PyRosetta protocol."""

        if clients_indices is None:
            return 0
        else:
            _protocols_index = len(clients_indices) - len(protocols)
            return clients_indices[_protocols_index]

    def _get_resource(
        self, resources: TaskChainResources, protocols: Sized
    ) -> TaskResource:
        """Return the resource for the current PyRosetta protocol."""

        if resources is None:
            return None
        else:
            _protocols_index = len(resources) - len(protocols)
            return resources[_protocols_index]

    def _get_priority(
        self, priorities: TaskChainPriorities, protocols: Sized
    ) -> Optional[int]:
        """Return the priority for the current PyRosetta protocol."""

        if priorities is None:
            return None
        else:
            _protocols_index = len(priorities) - len(protocols)
            return priorities[_protocols_index]

    def _get_retry(
        self, retries: TaskChainRetries, protocols: Sized
    ) -> Optional[int]:
        """Return the number of task retries for the current PyRosetta protocol."""

        if retries is None:
            return None
        elif isinstance(retries, int):
            return retries
        else:
            _protocols_index = len(retries) - len(protocols)
            return retries[_protocols_index]

    def _parse_resources(self, resources: Any) -> Any:
        """Parse the `resources` keyword argument value of the `PyRosettaCluster.distribute` method."""

        _dask_version_threshold = (2, 1, 0)
        if __dask_version__ < _dask_version_threshold:
            if resources is not None:
                _dask_version_str = _version_tuple_to_str(__dask_version__)
                _dask_version_threshold_str = _version_tuple_to_str(_dask_version_threshold)
                logging.warning(
                    "Use of the `resources` keyword argument is not supported for 'dask' and 'distributed' "
                    f"package versions <{_dask_version_threshold_str}\nCurrent Dask version: {_dask_version_str}\n"
                    "Please set `PyRosettaCluster().distribute(resources=None)`, or upgrade the 'dask' and 'distributed' "
                    f"package versions to >={_dask_version_threshold_str} to silence this warning. "
                    "Automatically disabling resource constraints..."
                )
            return None
        else:
            return resources

    def _parse_priorities(self, priorities: Any) -> Any:
        """Parse the `priorities` keyword argument value of the `PyRosettaCluster.distribute` method."""

        _dask_version_threshold = (1, 21, 0)
        if __dask_version__ < _dask_version_threshold:
            if priorities is not None:
                _dask_version_str = _version_tuple_to_str(__dask_version__)
                _dask_version_threshold_str = _version_tuple_to_str(_dask_version_threshold)
                logging.warning(
                    "Use of the `priorities` keyword argument is not supported for 'dask' and 'distributed' "
                    f"package versions <{_dask_version_threshold_str}\nCurrent Dask version: {_dask_version_str}\n"
                    "Please set `PyRosettaCluster().distribute(priorities=None)`, or upgrade the 'dask' and 'distributed' "
                    f"package versions to >={_dask_version_threshold_str} to silence this warning. "
                    "Automatically disabling priorities..."
                )
            return None
        else:
            return priorities

    def _parse_retries(self, retries: Any) -> Any:
        """Parse the `retries` keyword argument value of the `PyRosettaCluster.distribute` method."""

        _dask_version_threshold = (1, 20, 0)
        if __dask_version__ < _dask_version_threshold:
            if retries is not None:
                _dask_version_str = _version_tuple_to_str(__dask_version__)
                _dask_version_threshold_str = _version_tuple_to_str(_dask_version_threshold)
                logging.warning(
                    "Use of the `retries` keyword argument is not supported for 'dask' and 'distributed' "
                    f"package versions <{_dask_version_threshold_str}\nCurrent Dask version: {_dask_version_str}\n"
                    "Please set `PyRosettaCluster().distribute(retries=None)`, or upgrade the 'dask' and 'distributed' "
                    f"package versions to >={_dask_version_threshold_str} to silence this warning. "
                    "Automatically disabling task retries..."
                )
            return None
        else:
            return retries

    def _setup_kwargs(
        self,
        kwargs: Dict[str, Any],
        clients_indices: TaskChainClientIndices,
        resources: TaskChainResources,
        priorities: TaskChainPriorities,
        retries: TaskChainRetries,
    ) -> Tuple[
        bytes,
        Dict[str, Any],
        PyRosettaProtocol,
        int,
        TaskResource,
        Optional[int],
        Optional[int],
    ]:
        """Setup the keyword arguments for the subsequent user-defined task dictionary."""

        clients_index = self._get_clients_index(clients_indices, kwargs[self.protocols_key])
        resource = self._get_resource(resources, kwargs[self.protocols_key])
        priority = self._get_priority(priorities, kwargs[self.protocols_key])
        retry = self._get_retry(retries, kwargs[self.protocols_key])
        _protocols, protocol, seed = self._get_task_state(kwargs[self.protocols_key])
        kwargs[self.protocols_key] = _protocols
        kwargs = self._setup_seed(kwargs, seed)
        pyrosetta_init_kwargs = self._setup_pyrosetta_init_kwargs(kwargs)
        compressed_kwargs = self.serializer.compress_kwargs(kwargs)

        return (
            compressed_kwargs,
            pyrosetta_init_kwargs,
            protocol,
            clients_index,
            resource,
            priority,
            retry,
        )

    def _setup_seed(self, kwargs: Dict[str, Any], seed: Optional[str]) -> Dict[str, Any]:
        """
        Update the value of the "options" or "extra_options" key of the user-defined task dictionary with the
        `-run:jran` Rosetta command-line option.
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
        self,
        args: Tuple[Any, ...],
        protocols: Any,
        clients_indices: Any,
        resources: Any,
        priorities: Any,
        retries: Any,
    ) -> Tuple[
        PyRosettaProtocols,
        PyRosettaProtocol,
        Optional[str],
        int,
        TaskResource,
        Optional[int],
        Optional[int],
    ]:
        """Parse, validate, and setup the user-defined PyRosetta protocol(s)."""

        _protocols = _parse_protocols(args) + _parse_protocols(protocols)
        _clients_dict_keys = list(self.clients_dict.keys())
        _validate_clients_indices(clients_indices, _protocols, _clients_dict_keys)
        _validate_resources(resources, _protocols)
        _validate_priorities(priorities, _protocols)
        _validate_retries(retries, _protocols)
        _clients_index = self._get_clients_index(clients_indices, _protocols)
        _resource = self._get_resource(resources, _protocols)
        _priority = self._get_priority(priorities, _protocols)
        _retry = self._get_retry(retries, _protocols)
        _protocols, _protocol, _seed = self._get_task_state(
            _validate_protocols_seeds_decoy_ids(_protocols, self.seeds, self.decoy_ids)
        )

        return (
            _protocols,
            _protocol,
            _seed,
            _clients_index,
            _resource,
            _priority,
            _retry,
        )


def capture_task_metadata(func: CallableType) -> CallableType:
    """Capture metadata into the user-defined task dictionary."""

    @wraps(func)
    def wrapper(
        protocol_name: str,
        protocol: PyRosettaProtocol,
        packed_pose: PackedPose,
        datetime_format: str,
        norm_task_options: bool,
        ignore_errors: bool,
        protocols_key: str,
        decoy_ids: List[int],
        serializer: Serialization,
        kwargs: Dict[str, Any],
    ) -> Any:
        """Wrapper function to the `capture_task_metadata` decorator."""

        kwargs["PyRosettaCluster_protocol_name"] = protocol_name
        if not "PyRosettaCluster_protocols" in kwargs:
            kwargs["PyRosettaCluster_protocols"] = []
        kwargs["PyRosettaCluster_protocols"].append(protocol_name)
        protocol_number = (
            len(kwargs["PyRosettaCluster_protocols"]) - 1
        )  # Pythonic indices
        kwargs["PyRosettaCluster_protocol_number"] = protocol_number
        if not "PyRosettaCluster_datetime_start" in kwargs:
            kwargs["PyRosettaCluster_datetime_start"] = datetime.now().strftime(
                datetime_format
            )
        seed = pyrosetta.rosetta.numeric.random.rg().get_seed()
        if not "PyRosettaCluster_seeds" in kwargs:
            kwargs["PyRosettaCluster_seeds"] = []
        kwargs["PyRosettaCluster_seeds"].append((protocol_name, str(seed)))
        kwargs["PyRosettaCluster_seed"] = seed
        if norm_task_options and protocol_number == 0:
            options = _get_norm_task_options(ignore_errors)
            if "extra_options" in kwargs["PyRosettaCluster_task"]:
                kwargs["PyRosettaCluster_task"]["extra_options"] = options
                if "options" in kwargs["PyRosettaCluster_task"]:
                    kwargs["PyRosettaCluster_task"]["options"] = ""
            elif "options" in kwargs["PyRosettaCluster_task"]:
                kwargs["PyRosettaCluster_task"]["options"] = options
            else:
                kwargs["PyRosettaCluster_task"]["extra_options"] = options

        return func(
            protocol_name,
            protocol,
            packed_pose,
            datetime_format,
            norm_task_options,
            ignore_errors,
            protocols_key,
            decoy_ids,
            serializer,
            kwargs,
        )

    return cast(CallableType, wrapper)
