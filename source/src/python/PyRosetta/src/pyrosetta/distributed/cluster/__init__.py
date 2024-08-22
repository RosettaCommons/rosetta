# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

import warnings

from contextlib import suppress
from pyrosetta.distributed.cluster.core import PyRosettaCluster
from pyrosetta.distributed.cluster.toolkit import (
    Serialization,
    get_instance_kwargs,
    get_protocols,
    get_protocols_list_of_str,
    get_scores_dict,
    get_yml,
    produce,
    recreate_environment,
    reproduce,
    reserve_scores,
    run,
    update_scores,
    _print_conda_warnings,
)
from typing import List


__all__: List[str] = [
    "PyRosettaCluster",
    "Serialization",
    "get_instance_kwargs",
    "get_protocols",
    "get_protocols_list_of_str",
    "get_scores_dict",
    "get_yml",
    "produce",
    "recreate_environment",
    "reproduce",
    "reserve_scores",
    "run",
    "update_scores",
]
__version__: str = "1.2.1"

_print_conda_warnings()

with warnings.catch_warnings() and suppress(NameError):
    warnings.simplefilter("ignore", category=UserWarning)
    # Catch warning inside ipython interpreter:
    #     UserWarning: `IPython.core.IPCompleter.limit_to__all__` configuration
    #     value has been deprecated since IPython 5.0, will be made to have no
    #     effects and then removed in future version of IPython.
    # Suppress exception outside ipython interpreter:
    #     NameError: name 'get_ipython' is not defined
    get_ipython().Completer.limit_to__all__ = True
