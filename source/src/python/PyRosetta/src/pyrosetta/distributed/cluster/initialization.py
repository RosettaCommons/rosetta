# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

try:
    import numpy
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.initialization' requires the "
        + "third-party packages 'numpy' and 'toolz' as dependencies!\n"
        + "Please install the packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/numpy/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import inspect
import logging
import os
import pyrosetta
import pyrosetta.distributed

from pyrosetta.exceptions import PyRosettaIsNotInitializedError
from pyrosetta.rosetta.basic import was_init_called
from typing import (
    AbstractSet,
    Dict,
    List,
    NoReturn,
    Optional,
    Union,
)


def _get_pyrosetta_init_args() -> List[str]:
    """Return a `list` object representing `pyrosetta.init` parameters."""
    return inspect.getfullargspec(pyrosetta.init).args


def _get_residue_type_set_name3() -> Union[AbstractSet[str], NoReturn]:
    """
    Return a `set` of `str` of 3-letter names of residues in the PyRosetta ResidueTypeSet database.
    """

    if was_init_called():
        _pose = pyrosetta.Pose()
        _res_set = _pose.conformation().modifiable_residue_type_set_for_conf()
        _brt = _res_set.base_residue_types()
        _name3_list = [_brt.pop().name3() for _ in range(_brt.capacity())]
    else:
        raise PyRosettaIsNotInitializedError("PyRosetta is not initialized.")

    return set(_name3_list)


def _maybe_init_client() -> Optional[NoReturn]:
    """
    Initialize PyRosetta if it has not been initialized, otherwise confirm that PyRosetta was initialized with
    a constant seed.
    """

    if was_init_called():
        err_msg = (
            "Critical error! It appears that PyRosetta was already initialized without "
            "the '-run:constant_seed 1' option! Therefore, any work done "
            "before instantiating `PyRosettaCluster` cannot be reproduced! "
            "`PyRosettaCluster` is a tool for reproducible macromolecular modeling and design. "
            "If you're running `PyRosettaCluster` in a Jupyter Notebook or JupyterLab, "
            "please restart the kernel and initialize PyRosetta with the Rosetta command-line option "
            "'-run:constant_seed 1' before preparing any input `Pose` or `PackedPose` objects; "
            "i.e., run `pyrosetta.init('-run:constant_seed 1')`. If you are running from a Python script, "
            "please initialize PyRosetta with the Rosetta command-line option '-run:constant_seed 1'. "
            "The following actual PyRosetta RNG seed was retrieved compared to the desired constant seed:"
        )
        numpy.testing.assert_equal(
            pyrosetta.rosetta.numeric.random.rg().get_seed(),
            1111111,
            err_msg=err_msg,
            verbose=True,
        )
    else:
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-mute all",
        )


@toolz.functoolz.curry
def _maybe_relativize(
    option_name: str, value: str, start: str, ignore_errors: bool
) -> Union[str, NoReturn]:
    """Relativize a `str` object if it exists as a filesystem path."""

    try:
        expanded = os.path.expandvars(os.path.expanduser(value))
        maybe_path = os.path.normpath(
            expanded if os.path.isabs(expanded) else os.path.join(start, expanded)
        )
        return os.path.relpath(maybe_path, start=start) if os.path.lexists(maybe_path) else value
    except Exception as ex: # May be malformed object or cross-drive path on Windows
        _msg = (
            "`PyRosettaCluster` with the `norm_task_options` instance attribute set to `True` cannot "
            f"relativize a path in the task's Rosetta command-line option '-{option_name}' with value: '{value}'"
        )
        if ignore_errors:
            logging.error(
                (
                    "{0}: {1}. {2}. Please consider setting a relative path in the task's Rosetta "
                    "command-line option value. Ignoring error because the `ignore_errors` instance "
                    "attribute is set to `True`!"
                ).format(type(ex).__name__, ex, _msg)
            )
            return value
        else:
            raise ValueError(
                (
                    "{0}. {1}. Please update the task's Rosetta command-line option value, "
                    "set the `norm_task_options` keyword argument value to `False`, or set the "
                    "`ignore_errors` keyword argument value to `True` in order to continue."
                ).format(ex, _msg)
            )


def _get_norm_task_options(ignore_errors: bool) -> Dict[str, str]:
    """Get normalized Rosetta command-line options for a task dictionary."""

    options_dict: Dict[str, List[str]] = toolz.dicttoolz.keyfilter(
        lambda k: k not in ("in:path:database", "run:constant_seed", "run:jran"),
        pyrosetta.get_init_options(compressed=False, as_dict=True),
    )
    relativize = _maybe_relativize(start=os.getcwd(), ignore_errors=ignore_errors)

    msgs: List[str] = []
    options: Dict[str, str] = {}
    for option_name, values in options_dict.items():
        rel_values = []
        for value in values:
            rel_value = relativize(option_name, value)
            if rel_value != value:
                msgs.append(f"'{value}' -> '{rel_value}'")
            rel_values.append(rel_value)
        options[option_name] = " ".join(rel_values)

    if msgs:
        logging.info(
            "`PyRosettaCluster` with the `norm_task_options` instance attribute set to `True` "
            + "is normalizing the following values in the task's Rosetta command-line options:\n"
            + "\n".join(msgs)
        )

    return options
