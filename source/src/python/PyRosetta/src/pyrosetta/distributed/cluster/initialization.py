# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import numpy as np
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.initialization' requires the "
        + "third-party package 'numpy' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/numpy/\n"
    )
    raise

import inspect
import pyrosetta
import pyrosetta.distributed

from typing import (
    AbstractSet,
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
