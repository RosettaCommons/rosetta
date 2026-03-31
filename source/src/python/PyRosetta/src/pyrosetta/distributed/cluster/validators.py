# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

try:
    import dask
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.validators' requires the "
        + "third-party package 'dask' as a dependency!\n"
        + "Please install this package into your virtual environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
    )
    raise

import json
import logging
import os

from typing import (
    AbstractSet,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Optional,
    Union,
)

OPTION_KEYS: AbstractSet[str] = {"options", "extra_options"}
DISALLOWED_RUN_OPTIONS: AbstractSet[str] = set(
    f"{prefix}{option}"
    for prefix in ("-run::", "-run:", "-")
    for option in ("constant_seed", "jran", "use_time_as_seed", "rng_seed_device", "seed_offset", "rng")
)
DISALLOWED_RUN_OPTIONS_NO_PREFIX: AbstractSet[str] = set(
    map(lambda option: option.lstrip("-"), DISALLOWED_RUN_OPTIONS)
)
PYROSETTACLUSTER_KEY_PREFIX: str = "PyRosettaCluster_"


def _validate_clients_indices(
        clients_indices: Any, _protocols: List[Callable[..., Any]], _clients_dict_keys: List[int],
    ) -> None:
    """Validate the `clients_indices` keyword argument value for the `PyRosettaCluster.distribute` method."""

    if clients_indices is not None:
        if not isinstance(clients_indices, (tuple, list)):
            raise ValueError(
                "The `clients_indices` keyword argument value must be of type `list` or `tuple`.\n"
                + f"Received: {type(clients_indices)}\n"
            )
        for i in clients_indices:
            if not isinstance(i, int):
                raise ValueError(
                    "Each item of the `clients_indices` keyword argument value must be of type `int`.\n"
                    + f"Received: {type(i)}\n"
                )
        if len(clients_indices) != len(_protocols):
            raise ValueError(
                "The `clients_indices` keyword argument value must have the same length as the number of user-defined PyRosetta protocols!\n"
                + f"Received `PyRosettaCluster.distribute(protocols=...)`: {_protocols}\n"
                + f"Received `PyRosettaCluster.distribute(clients_indices=...)`: {clients_indices}\n"
            )
        if not all(x in _clients_dict_keys for x in set(clients_indices)):
            raise ValueError(
                "Each item of the `clients_indices` keyword argument value must correspond to an index passed to the `PyRosettaCluster(clients=...)` keyword argument value.\n"
                + f"Available clients indices based on the `PyRosettaCluster(clients=...)` instance attribute: {_clients_dict_keys}\n"
                + f"Received `PyRosettaCluster.distribute(clients_indices=...)`: {clients_indices}\n"
            )
        for i in _clients_dict_keys:
            if clients_indices.count(i) == 0:
                logging.warning(
                    f"The `Client` object at index {i} of the `PyRosettaCluster(clients=...)` instance attribute will not be used in the simulation "
                    + "because it was not added as an index in the `PyRosettaCluster().distribute(clients_indices=...)` keyword argument value!\n"
                    + f"Available clients indices based on the `PyRosettaCluster(clients=...)` instance attribute: {_clients_dict_keys}\n"
                    + f"Received `PyRosettaCluster().distribute(clients_indices=...)`: {clients_indices}\n"
                    + "Continuing with the simulation.\n"
                )


def _validate_resources(resources: Any, _protocols: List[Callable[..., Any]]) -> None:
    """Validate the `resources` keyword argument value for the `PyRosettaCluster.distribute` method."""

    if resources is not None:
        if not isinstance(resources, (tuple, list)):
            raise ValueError(
                "The `resources` keyword argument value must be of type `list` or `tuple`. "
                + f"Received: {type(resources)}"
            )
        for obj in resources:
            if not isinstance(obj, (dict, type(None))):
                raise ValueError(
                    "Each item of the `resources` keyword argument value must be of type `dict` or `NoneType`. "
                    + f"Received: {type(obj)}"
                )
        if len(resources) != len(_protocols):
            raise ValueError(
                "The `resources` keyword argument value must have the same length as the number of user-defined PyRosetta protocols!\n"
                + f"Received `PyRosettaCluster().distribute(protocols=...)`: {_protocols}\n"
                + f"Received `PyRosettaCluster().distribute(resources=...)`: {resources}\n"
            )


def _validate_priorities(priorities: Any, _protocols: List[Callable[..., Any]]) -> None:
    """Validate the `priorities` keyword argument value for the `PyRosettaCluster.distribute` method."""

    if priorities is not None:
        if not isinstance(priorities, (tuple, list)):
            raise ValueError(
                "The `priorities` keyword argument value must be of type `list` or `tuple`. "
                + f"Received: {type(priorities)}"
            )
        for obj in priorities:
            if not isinstance(obj, int):
                raise ValueError(
                    "Each item of the `priorities` keyword argument value must be of type `int`. "
                    + f"Received: {type(obj)}"
                )
        if len(priorities) != len(_protocols):
            raise ValueError(
                "The `priorities` keyword argument value must have the same length as the number of user-defined PyRosetta protocols!\n"
                + f"Received `PyRosettaCluster.distribute(protocols=...)`: {_protocols}\n"
                + f"Received `PyRosettaCluster.distribute(priorities=...)`: {priorities}\n"
            )


def _validate_retries(retries: Any, _protocols: List[Callable[..., Any]]) -> None:
    """Validate the `retries` keyword argument value for the `PyRosettaCluster.distribute` method."""

    if retries is not None:
        if isinstance(retries, int):
            if retries < 0:
                raise ValueError(
                    "If the `retries` keyword argument value is of type `int`, it must be greater than or equal to 0.\n"
                    + f"Received: {retries}\n"
                )
        elif isinstance(retries, (tuple, list)):
            for obj in retries:
                if not isinstance(obj, int):
                    raise ValueError(
                        "Each item of the `retries` keyword argument value must be of type `int`. "
                        + f"Received: {type(obj)}"
                    )
                if obj < 0:
                    raise ValueError(
                        "Each item of the `retries` keyword argument value must be greater than or equal to 0. "
                        + f"Received: {obj}"
                    )
            if len(retries) != len(_protocols):
                raise ValueError(
                    "The `retries` keyword argument value must have the same length as the number of user-defined PyRosetta protocols!\n"
                    + f"Received `PyRosettaCluster.distribute(protocols=...)`: {_protocols}\n"
                    + f"Received `PyRosettaCluster.distribute(retries=...)`: {retries}\n"
                )
        else:
            raise ValueError(
                "The `retries` keyword argument value must be of type `int`, `list`, or `tuple`. "
                + f"Received: {type(retries)}\n"
            )


def _validate_scorefile_name(self, attribute: str, value: Any) -> None:
    """Validate the `scorefile_name` keyword argument value of `PyRosettaCluster`."""

    if not value.endswith(".json"):
        raise ValueError(f"The '{attribute}' keyword argument value must end in '.json'. Received: '{value}'")


def _validate_output_init_file(self, attribute: str, value: Any) -> None:
    """Validate the `output_init_file` keyword argument value of `PyRosettaCluster`."""

    if value != "" and not value.endswith(".init"):
        raise ValueError(f"The '{attribute}' keyword argument value must end in '.init'. Received: '{value}'")


def _validate_logging_address(self, attribute: str, value: Any) -> None:
    """Validate the `logging_address` keyword argument value of `PyRosettaCluster`."""

    if not isinstance(value, str):
        raise ValueError(f"The `{attribute}` keyword argument value must be of type `str`. Received: {type(value)}")
    if value.count(":") != 1:
        raise ValueError(f"The `{attribute}` keyword argument value must contain one colon. Received: '{value}'")
    _host, _port = tuple(s.strip() for s in value.split(":"))
    if not _host:
        raise ValueError(
            f"The `{attribute}` keyword argument value must contain a value before the colon representing the host. "
            + f"Received: '{value}'"
        )
    if not _port.isdigit():
        raise ValueError(
            f"The `{attribute}` keyword argument value must contain a digit after the colon representing the port. "
            + f"Received: '{value}'"
        )


def _validate_dirs(self, attribute: str, value: Any) -> None:
    """Validate and make the output, logging, and decoy directories."""

    output_dir = os.path.abspath(self.output_path)
    decoy_dir = os.path.join(output_dir, self.decoy_dir_name)
    logs_dir = os.path.join(output_dir, self.logs_dir_name)
    for d in [output_dir, decoy_dir, logs_dir]:
        os.makedirs(d, exist_ok=True)
        if not os.path.isdir(d):
            raise ValueError(f"`{d}` directory could not be created.")


def _validate_dir(self, attribute: str, value: str) -> None:
    """Validate the scratch directory."""

    if value != "":
        os.makedirs(value, exist_ok=True)
        if not os.path.isdir(value):
            raise ValueError(f"The `{attribute}` directory {value} could not be created.")


def _validate_int(self, attribute: str, value: int) -> None:
    """Validate that integers are greater than or equal to 1."""

    if value < 1:
        raise ValueError(
            f"The `{attribute}` must be a positive integer greater than or equal to 1."
        )
    

def _validate_min_len(self, attribute: str, value: Optional[List[Any]]) -> None:
    """Optionally validate that iterables have at least one object."""

    if value is not None and len(value) < 1:
        raise ValueError(
            f"The `{attribute}` must have at least one item if not `None`."
        )


def _validate_max_task_replicas(self, attribute: str, value: Optional[int]) -> None:
    """
    Validate that the value is `None` or integers are greater than or equal to 0, and that Dask's Active Memory
    Manager (AMM) policy is disabled.
    """

    if not (value is None or (isinstance(value, int) and value >= 0)):
        raise ValueError(f"`{attribute}` must be `None` or a positive integer greater than or equal to 0. Received: {value}")
    if isinstance(value, int) and value > 0:
        amm_start = dask.config.get("distributed.scheduler.active-memory-manager.start")
        amm_policies = dask.config.get("distributed.scheduler.active-memory-manager.policies") or []
        reduce_replicas_enabled = any(
            policy.get("class") == "distributed.active_memory_manager.ReduceReplicas"
            for policy in amm_policies
        )
        if amm_start and reduce_replicas_enabled:
            raise ValueError(
                "To use task replicas, please (1) disable Dask's `ReduceReplicas` policy, or (2) disable the Active Memory Manager (AMM) entirely. "
                + "For (1), run the following (before instantiating `PyRosettaCluster` or any `distributed.Client` objects):\n"
                + "    dask.config.set({'distributed.scheduler.active-memory-manager.policies': [...]}\n"
                + "Ensure the `[...]` list does not contain the entry: `{'class': 'distributed.active_memory_manager.ReduceReplicas'}`"
                + "For (2), run the following (before instantiating `PyRosettaCluster` or any `distributed.Client` objects):\n"
                + "    dask.config.set({'distributed.scheduler.active-memory-manager.start': False})\n"
                + "For more information, see https://distributed.dask.org/en/stable/active_memory_manager.html#reducereplicas "
                + "and https://docs.dask.org/en/stable/configuration.html"
            )


def _validate_float(
    self, attribute: str, value: Union[float, int]
) -> None:
    """Validate that `float` or `int` objects are greater than or equal to 0."""

    msg = (
        f"The `{attribute}` keyword argument value must be a positive `float` or `int` value greater than "
        f"or equal to 0. Received: {value}"
    )
    try:
        float(value)
    except:
        raise ValueError(msg)
    if value < 0.0:
        raise ValueError(msg)


def _validate_protocols_seeds_decoy_ids(
    protocols: List[Union[Callable[..., Any], Iterable[Any]]],
    seeds: List[str],
    decoy_ids: List[int],
) -> List[Union[Callable[..., Any], Iterable[Any]]]:
    """
    Validate that the user-defined PyRosetta protocols, and the `seeds` and `decoy_ids` keyword argument values
    of `PyRosettaCluster` have the same size.
    """

    if len(protocols) < 1:
        raise RuntimeError(
            "The user-defined PyRosetta protocols must contain at least one "
            + "callable or generator."
        )
    if seeds:
        assert len(protocols) == len(
            seeds
        ), "Number of input seeds must match number of input protocols!"
    if decoy_ids:
        assert len(protocols) == len(
            decoy_ids
        ), "Number of input decoy ids must match number of input protocols!"
    if seeds and decoy_ids:
        assert len(seeds) == len(
            decoy_ids
        ), "Number of input decoy ids must match number of input seeds!"

    return protocols


def _validate_residue_type_sets(
    _target_residue_type_set: Optional[AbstractSet[str]] = None,
    _client_residue_type_set: Optional[AbstractSet[str]] = None,
) -> None:
    """Validate that the ResidueType set in the `billiard` subprocess equals that in the head node process."""

    if _target_residue_type_set != _client_residue_type_set:
        _msg = (
            "The spawned subprocess ResidueType set does not equal the head node process ResidueType set! "
            "Please ensure that the spawned subprocess and the head node process have initialized PyRosetta "
            "with identical '-extra_res_fa', '-extra_res_cen', '-extra_patch_fa' and '-extra_patch_cen' "
            "Rosetta command-line options.\n"
            f"Spawned subprocess unique ResidueTypes: {_target_residue_type_set.difference(_client_residue_type_set)}\n"
            f"Head node process unique ResidueTypes: {_client_residue_type_set.difference(_target_residue_type_set)}\n"
        )
        logging.error(_msg)
        raise AssertionError(_msg)


def _validate_task(task: Dict[Any, Any]) -> None:
    """Validate that a task does not contain disallowed or reserved keys/values."""

    _msg = (
        "Disallowed Rosetta command-line option '{option}' in the value of the '{key}' key of task: {task}\n"
        "`PyRosettaCluster` handles seeding automatically. Please remove this Rosetta command-line option to continue."
    )
    for k, v in task.items():
        if not isinstance(k, str):
            raise ValueError(
                f"User-defined task dictionary key must be an instance of `str`. Received {type(k)} in task: {task}"
            )
        elif k in OPTION_KEYS:
            if isinstance(v, dict):
                for _option in v.keys():
                    if any(_option in x for x in (DISALLOWED_RUN_OPTIONS, DISALLOWED_RUN_OPTIONS_NO_PREFIX)):
                        raise ValueError(_msg.format(option=_option, key=k, task=task))
            elif isinstance(v, str):
                for _disallowed_option in DISALLOWED_RUN_OPTIONS:
                    if _disallowed_option in v:
                        raise ValueError(_msg.format(option=_disallowed_option, key=k, task=task))
            else:
                raise ValueError(
                    f"The value of the '{k}' key must be an instance of `dict` or `str`. Received {type(v)} in task: {task}"
                )
        elif k.startswith(PYROSETTACLUSTER_KEY_PREFIX):
            raise ValueError(
                f"Disallowed user-defined task dictionary key '{k}' of task: {task}\n"
                + f"Task keys starting with '{PYROSETTACLUSTER_KEY_PREFIX}' are reserved for PyRosettaCluster."
            )


def _validate_tasks(self, attribute: str, value: List[Dict[Any, Any]]) -> None:
    """Validate that tasks do not contain disallowed or reserved keys/values."""

    for task in value:
        _validate_task(task)
        if not _is_json_roundtrip_equal(task):
            raise ValueError(
                "An input user-defined task dictionary is not JSON-serializable! Please reformat "
                + f"the input user-defined task dictionary for the PyRosettaCluster simulation:\n{task}"
            )


def _is_json_roundtrip_equal(obj: Any) -> bool:
    """Test if an object is equal after roundtrip JSON-serialization."""

    try:
        return obj == json.loads(json.dumps(obj))
    except (TypeError, json.JSONDecodeError, ValueError):
        return False
    except Exception as ex:
        raise RuntimeError("Unexpected JSON-serialization error during roundtrip validation.") from ex
