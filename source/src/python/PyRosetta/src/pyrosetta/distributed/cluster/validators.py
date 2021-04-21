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


def _validate_dirs(self, attribute: str, value: Any) -> None:
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
    _master_residue_type_set: Optional[AbstractSet[str]] = None,
) -> Optional[NoReturn]:
    """
    Validate that the compute instance (distributed worker) ResidueType set equals
    master instance (local host) ResidueType set.
    """

    if _target_residue_type_set != _master_residue_type_set:
        _msg = (
            "Error! The compute instance (distributed worker) ResidueType set does not equal "
            + "the master instance (local host) ResidueType set! Please ensure that the "
            + "compute instance (distributed worker) and master instance (local host) have "
            + "initialized PyRosetta with identical '-extra_res_fa' and '-extra_res_cen' options."
        )
        logging.error(_msg)
        raise AssertionError(_msg)
