# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

import os
import shutil
import sys

from pyrosetta.utility import get_package_version
from typing import Dict, Generic, List, Tuple, TypeVar


__dask_version__: Tuple[int, int, int] = get_package_version("dask")
__dask_jobqueue_version__: Tuple[int, int, int] = get_package_version("dask-jobqueue")

G = TypeVar("G")


class EnvironmentConfig(Generic[G]):
    _env_var: str = "PYROSETTACLUSTER_ENVIRONMENT_MANAGER"
    _environment_cmds: Dict[str, str] = {
        "pixi": "pixi workspace export conda-environment",
        "uv": "uv export --format requirements-txt --frozen",
        "mamba": "mamba env export --prefix {0}".format(sys.prefix),
        "conda": "conda env export --prefix {0}".format(sys.prefix),
    }

    def __init__(self) -> None:
        _env_var_manager = os.environ.get(EnvironmentConfig._env_var, None)
        if _env_var_manager:
            self.environment_manager = _env_var_manager
            print(
                "Configuring environment manager for PyRosettaCluster from operating system "
                + f"environment variable: {EnvironmentConfig._env_var}={self.environment_manager}"
            )
            if self.environment_manager not in EnvironmentConfig._environment_cmds.keys():
                raise ValueError(
                    "The '{0}' environment variable must be in: {1}. Received: '{2}'. ".format(
                        EnvironmentConfig._env_var,
                        list(EnvironmentConfig._environment_cmds.keys()),
                        self.environment_manager,
                    )
                )
        else:
            for _manager in EnvironmentConfig._environment_cmds.keys():
                if shutil.which(_manager):
                    self.environment_manager = _manager
                    print(f"Configuring environment manager for PyRosettaCluster: '{_manager}'")
                    break
            else:
                self.environment_manager = "conda"
                print(
                    f"Warning: could not configure an environment manager for PyRosettaCluster. "
                    + "Please ensure that either of 'pixi', 'uv', 'mamba', or 'conda' is installed. "
                    + "Using 'conda' as the default environment manager."
                )
        self.environment_cmd = self._environment_cmds[self.environment_manager]


_env_config = EnvironmentConfig()
environment_manager: str = _env_config.environment_manager
environment_cmd: str = _env_config.environment_cmd

source_domains: List[str] = [
    "conda.graylab.jhu.edu",
    "west.rosettacommons.org",
    "conda.rosettacommons.org",
]  # Conda channels and/or source domains (containing PyRosetta usernames/passwords) to be stripped from YML file strings.
