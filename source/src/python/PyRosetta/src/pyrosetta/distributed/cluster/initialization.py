# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import numpy as np
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

from typing import (
    AbstractSet,
    Dict,
    List,
    NoReturn,
    Optional,
    Union,
)


def _get_pyrosetta_init_args() -> List[str]:
    """
    Return a `list` object representing `pyrosetta.init` arguments.
    """
    return inspect.getfullargspec(pyrosetta.init).args


def _get_residue_type_set_name3() -> Union[AbstractSet[str], NoReturn]:
    """
    Return a `set` of `str` of 3-letter names of residues in the PyRosetta
    ResidueTypeSet database.
    """

    if pyrosetta.rosetta.basic.was_init_called():
        _pose = pyrosetta.Pose()
        _res_set = _pose.conformation().modifiable_residue_type_set_for_conf()
        _brt = _res_set.base_residue_types()
        _name3_list = [_brt.pop().name3() for _ in range(_brt.capacity())]
    else:
        raise RuntimeError("PyRosetta is not initialized.")

    return set(_name3_list)


def _maybe_init_client() -> Optional[NoReturn]:
    """
    Initialize PyRosetta if it has not been initialized, otherwise confirm that
    PyRosetta was initialized with a constant seed.
    """

    if pyrosetta.rosetta.basic.was_init_called():
        err_msg = (
            "Critical error! It appears that pyrosetta was already initialized without "
            "the '-run:constant_seed 1' option! Therefore, any work done "
            "before instantiating PyRosettaCluster cannot be reproduced! "
            "PyRosettaCluster is a tool for reproducible computational protein design. "
            "If you're running PyRosettaCluster in a Jupyter Notebook or Jupyter Lab, "
            "please restart the kernel and initialize pyrosetta with the command line option "
            "'-run:constant_seed 1' before preparing any input poses; "
            "i.e. pyrosetta.init('-run:constant_seed 1'). If you're running from a python script, "
            "please initialize pyrosetta with command line option '-run:constant_seed 1'. "
            "The following actual seed was retrieved compared to the desired constant seed:"
        )
        np.testing.assert_equal(
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
    """Relativize a `str` object if it exists as a path."""

    _log_msg = (
        "{0}: {1}. {2}. Please consider setting a relative path in the task's PyRosetta "
        "initialization option value. Ignoring error because `ignore_errors` is enabled!"
    )
    _err_msg = (
        "{0}. {1}. Please update the task's PyRosetta initialization option value, "
        "disable `norm_task_options`, or enable `ignore_errors` to continue."
    )

    try:
        expanded = os.path.expandvars(os.path.expanduser(value))
        maybe_path = os.path.normpath(
            expanded if os.path.isabs(expanded) else os.path.join(start, expanded)
        )
    except Exception as ex: # May be malformed object
        msg = (
            "PyRosettaCluster (with `norm_task_options` enabled) cannot construct a candidate path "
            f"in the task's PyRosetta initialization option '-{option_name}' with value: '{value}'"
        )
        if ignore_errors:
            logging.error(_log_msg.format(type(ex).__name__, ex, msg))
            return value
        else:
            raise ValueError(_err_msg.format(ex, msg))

    if os.path.lexists(maybe_path):
        try:
            return os.path.relpath(maybe_path, start=start)
        except Exception as ex: # May be cross-drive path on Windows
            msg = (
                "PyRosettaCluster (with `norm_task_options` enabled) cannot relativize a path in the "
                f"task's PyRosetta initialization option '-{option_name}' with value: '{value}'"
            )
            if ignore_errors:
                logging.error(_log_msg.format(type(ex).__name__, ex, msg))
                return value
            else:
                raise ValueError(_err_msg.format(ex, msg))
    else:
        return value


def _get_norm_task_options(ignore_errors: bool) -> Dict[str, str]:
    """Get normalized task PyRosetta initialization options."""

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
            "PyRosettaCluster (with `norm_task_options` enabled) is normalizing the "
            + "following values in the task's PyRosetta initialization options:\n"
            + "\n".join(msgs)
        )

    return options
