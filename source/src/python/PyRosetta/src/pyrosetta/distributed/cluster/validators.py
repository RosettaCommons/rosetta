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
    Dict,
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

    os.makedirs(value, exist_ok=True)
    if not os.path.isdir(value):
        raise ValueError(f"`scratch_dir` directory {value} could not be created.")


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
            "Error! The compute instance (distributed worker) ResidueType set does not equal "
            + "the client instance (local host) ResidueType set! Please ensure that the "
            + "compute instance (distributed worker) and client instance (local host) have "
            + "initialized PyRosetta with identical '-extra_res_fa' and '-extra_res_cen' options."
        )
        logging.error(_msg)
        raise AssertionError(_msg)
