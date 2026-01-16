# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

import logging
import os

from typing import (
    AbstractSet,
    Any,
    Callable,
    Iterable,
    List,
    NoReturn,
    Optional,
    Union,
)


def _validate_clients_indices(
        clients_indices: Any, _protocols: List[Callable[..., Any]], _clients_dict_keys: List[int],
    ) -> Optional[NoReturn]:
    """Validate the `clients_indices` keyword argument parameter."""
    if clients_indices is not None:
        if not isinstance(clients_indices, (tuple, list)):
            raise ValueError(
                "The `clients_indices` keyword argument parameter must be of type `list` or `tuple`.\n"
                + f"Received: {type(clients_indices)}\n"
            )
        for i in clients_indices:
            if not isinstance(i, int):
                raise ValueError(
                    "Each item of the `clients_indices` keyword argument parameter must be of type `int`.\n"
                    + f"Received: {type(i)}\n"
                )
        if len(clients_indices) != len(_protocols):
            raise ValueError(
                "The `clients_indices` keyword argument parameter must have the same length as the number of user-defined PyRosetta protocols!\n"
                + f"Received `PyRosettaCluster().distribute(protocols=...)`: {_protocols}\n"
                + f"Received `PyRosettaCluster().distribute(clients_indices=...)`: {clients_indices}\n"
            )
        if not all(x in _clients_dict_keys for x in set(clients_indices)):
            raise ValueError(
                "Each item of the `clients_indices` keyword argument parameter must correspond to an index passed to the `PyRosettaCluster(clients=...)` class attribute.\n"
                + f"Available clients indices based on the `PyRosettaCluster(clients=...)` class attribute: {_clients_dict_keys}\n"
                + f"Received `PyRosettaCluster().distribute(clients_indices=...)`: {clients_indices}\n"
            )
        for i in _clients_dict_keys:
            if clients_indices.count(i) == 0:
                logging.warning(
                    f"The `Client` object at index {i} of the `PyRosettaCluster(clients=...)` class attribute will not be used in the simulation "
                    + "because it was not added as an index in the `PyRosettaCluster().distribute(clients_indices=...)` keyword argument parameter!\n"
                    + f"Available clients indices based on the `PyRosettaCluster(clients=...)` class attribute: {_clients_dict_keys}\n"
                    + f"Received `PyRosettaCluster().distribute(clients_indices=...)`: {clients_indices}\n"
                    + "Continuing with the simulation.\n"
                )


def _validate_resources(resources: Any, _protocols: List[Callable[..., Any]]) -> Optional[NoReturn]:
    """Validate the `resources` keyword argument parameter."""
    if resources is not None:
        if not isinstance(resources, (tuple, list)):
            raise ValueError(
                "The `resources` keyword argument parameter must be of type `list` or `tuple`.\n"
                + f"Received: {type(resources)}\n"
            )
        for obj in resources:
            if not isinstance(obj, (dict, type(None))):
                raise ValueError(
                    "Each item of the `resources` keyword argument parameter must be of type `dict` or `NoneType`.\n"
                    + f"Received: {type(obj)}\n"
                )
        if len(resources) != len(_protocols):
            raise ValueError(
                "The `resources` keyword argument parameter must have the same length as the number of user-defined PyRosetta protocols!\n"
                + f"Received `PyRosettaCluster().distribute(protocols=...)`: {_protocols}\n"
                + f"Received `PyRosettaCluster().distribute(resources=...)`: {resources}\n"
            )


def _validate_priorities(priorities: Any, _protocols: List[Callable[..., Any]]) -> Optional[NoReturn]:
    """Validate the `priorities` keyword argument parameter."""
    if priorities is not None:
        if not isinstance(priorities, (tuple, list)):
            raise ValueError(
                "The `priorities` keyword argument parameter must be of type `list` or `tuple`.\n"
                + f"Received: {type(priorities)}\n"
            )
        for obj in priorities:
            if not isinstance(obj, int):
                raise ValueError(
                    "Each item of the `priorities` keyword argument parameter must be of type `int`.\n"
                    + f"Received: {type(obj)}\n"
                )
        if len(priorities) != len(_protocols):
            raise ValueError(
                "The `priorities` keyword argument parameter must have the same length as the number of user-defined PyRosetta protocols!\n"
                + f"Received `PyRosettaCluster().distribute(protocols=...)`: {_protocols}\n"
                + f"Received `PyRosettaCluster().distribute(priorities=...)`: {priorities}\n"
            )


def _validate_retries(retries: Any, _protocols: List[Callable[..., Any]]) -> Optional[NoReturn]:
    """Validate the `retries` keyword argument parameter."""
    if retries is not None:
        if isinstance(retries, int):
            if retries < 0:
                raise ValueError(
                    "If the `retries` keyword argument parameter is of type `int`, it must be greater than or equal to 0.\n"
                    + f"Received: {retries}\n"
                )
        elif isinstance(retries, (tuple, list)):
            for obj in retries:
                if not isinstance(obj, int):
                    raise ValueError(
                        "Each item of the `retries` keyword argument parameter must be of type `int`.\n"
                        + f"Received: {type(obj)}\n"
                    )
                if obj < 0:
                    raise ValueError(
                        "Each item of the `retries` keyword argument parameter must be greater than or equal to 0.\n"
                        + f"Received: {obj}\n"
                    )
            if len(retries) != len(_protocols):
                raise ValueError(
                    "The `retries` keyword argument parameter must have the same length as the number of user-defined PyRosetta protocols!\n"
                    + f"Received `PyRosettaCluster().distribute(protocols=...)`: {_protocols}\n"
                    + f"Received `PyRosettaCluster().distribute(retries=...)`: {retries}\n"
                )
        else:
            raise ValueError(
                "The `retries` keyword argument parameter must be of type `int`, `list`, or `tuple`.\n"
                + f"Received: {type(retries)}\n"
            )


def _validate_scorefile_name(self, attribute: str, value: Any) -> Optional[NoReturn]:
    if not value.endswith(".json"):
        raise ValueError(f"The '{attribute}' keyword argument parameter must end in '.json'.")


def _validate_output_init_file(self, attribute: str, value: Any) -> Optional[NoReturn]:
    if value != "" and not value.endswith(".init"):
        raise ValueError(f"The '{attribute}' keyword argument parameter must end in '.init'.")


def _validate_logging_address(self, attribute: str, value: Any) -> Optional[NoReturn]:
    if not isinstance(value, str):
        raise ValueError(f"`{attribute}` must be of type `str`. Received: '{type(value)}'")
    if value.count(":") != 1:
        raise ValueError(f"`{attribute}` must contain one colon. Received: '{value}'")
    _host, _port = tuple(s.strip() for s in value.split(":"))
    if not _host:
        raise ValueError(f"`{attribute}` must contain a value before the colon representing the host. Received: '{value}'")
    if not _port.isdigit():
        raise ValueError(f"`{attribute}` must contain a digit after the colon representing the port. Received: '{value}'")


def _validate_dirs(self, attribute: str, value: Any) -> Optional[NoReturn]:
    """Validate the output, logging, and decoy directories."""

    output_dir = os.path.abspath(self.output_path)
    decoy_dir = os.path.join(output_dir, self.decoy_dir_name)
    logs_dir = os.path.join(output_dir, self.logs_dir_name)
    for d in [output_dir, decoy_dir, logs_dir]:
        os.makedirs(d, exist_ok=True)
        if not os.path.isdir(d):
            raise ValueError(f"`{d}` directory could not be created.")


def _validate_dir(self, attribute: str, value: str) -> Optional[NoReturn]:
    """Validate the scratch directory."""
    if value != "":
        os.makedirs(value, exist_ok=True)
        if not os.path.isdir(value):
            raise ValueError(f"`{attribute}` directory {value} could not be created.")


def _validate_int(self, attribute: str, value: int) -> Optional[NoReturn]:
    """Validate that integers are greater than or equal to 1."""

    if value < 1:
        raise ValueError(
            f"`{attribute}` must be a positive integer greater than or equal to 1."
        )
    

def _validate_min_len(self, attribute: str, value: Optional[List[Any]]) -> Optional[NoReturn]:
    """Optionally validate that iterables have at least one object."""

    if value is not None and len(value) < 1:
        raise ValueError(
            f"`{attribute}` must have at least one item if not `None`."
        )


def _validate_float(
    self, attribute: str, value: Union[float, int]
) -> Optional[NoReturn]:
    """Validate that floats are greater than or equal to 0.0"""

    msg = f"`{attribute}` must be a positive `float` or `int` value greater than or equal to 0."
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
) -> Union[List[Union[Callable[..., Any], Iterable[Any]]], NoReturn]:
    """
    Validate that the user-provided PyRosetta protocols and PyRosettaCluster
    `seeds` and `decoy_ids` attributes have the same size.
    """
    if len(protocols) < 1:
        raise RuntimeError(
            "The user-provided PyRosetta protocols must contain at least one "
            + "function or generator."
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
) -> Optional[NoReturn]:
    """
    Validate that the compute instance (distributed worker) ResidueType set equals
    client instance (local host) ResidueType set.
    """

    if _target_residue_type_set != _client_residue_type_set:
        _msg = (
            "The compute instance (distributed worker) ResidueType set does not equal "
            + "the client instance (local host) ResidueType set! Please ensure that the "
            + "compute instance (distributed worker) and client instance (local host) have "
            + "initialized PyRosetta with identical '-extra_res_fa' and '-extra_res_cen' options.\n"
            + f"Compute instance unique ResidueTypes: {_target_residue_type_set.difference(_client_residue_type_set)}\n"
            + f"Client instance unique ResidueTypes: {_client_residue_type_set.difference(_target_residue_type_set)}\n"
        )
        logging.error(_msg)
        raise AssertionError(_msg)
