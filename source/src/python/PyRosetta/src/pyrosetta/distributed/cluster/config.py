# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    from distributed import get_worker
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.config' requires the "
        + "third-party package 'distributed' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/distributed/\n"
    )
    raise


import os
import shutil
import sys

from contextlib import contextmanager
from pyrosetta.utility import get_package_version
from typing import (
    Any,
    Dict,
    Generator,
    Generic,
    List,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
)


__dask_version__: Tuple[int, int, int] = get_package_version("dask")
__dask_jobqueue_version__: Tuple[int, int, int] = get_package_version("dask-jobqueue")

G = TypeVar("G")


class EnvironmentConfig(Generic[G]):
    _ENV_VAR: str = "PYROSETTACLUSTER_ENVIRONMENT_MANAGER"
    _ENV_MANAGERS: Tuple[str, ...] = ("pixi", "uv", "mamba", "conda")
    _ENV_EXPORT_CMDS: Dict[str, str] = {
        "pixi": "pixi workspace export conda-environment",
        "uv": "uv export --format requirements-txt --frozen",
        "mamba": "mamba env export --prefix {0}".format(sys.prefix),
        "conda": "conda env export --prefix {0}".format(sys.prefix),
    }
    _ENV_LIST_CMDS: Dict[str, Optional[str]] = {
        "pixi": None,
        "uv": None,
        "mamba": "mamba env list",
        "conda": "conda env list",
    }

    def __init__(self) -> None:
        _env_var_manager = os.environ.get(EnvironmentConfig._ENV_VAR, None)
        if _env_var_manager:
            self.environment_manager = _env_var_manager
            print(
                "Configuring environment manager for PyRosettaCluster from operating system "
                + f"environment variable: {EnvironmentConfig._ENV_VAR}={self.environment_manager}"
            )
            if self.environment_manager not in EnvironmentConfig._ENV_MANAGERS:
                raise ValueError(
                    "The '{0}' environment variable must be in: '{1}'. Received: '{2}'. ".format(
                        EnvironmentConfig._ENV_VAR,
                        EnvironmentConfig._ENV_MANAGERS,
                        self.environment_manager,
                    )
                )
        else:
            for _manager in EnvironmentConfig._ENV_MANAGERS:
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

    @property
    def env_export_cmd(self) -> str:
        return self._ENV_EXPORT_CMDS[self.environment_manager]

    @property
    def env_list_cmd(self) -> Optional[str]:
        return self._ENV_LIST_CMDS[self.environment_manager]

    def env_create_cmd(
        self, environment_name: str, raw_spec: str, tmp_dir: str
    ) -> Union[str, NoReturn]:
        if self.environment_manager in ("conda", "mamba", "pixi"):
            yml_file = os.path.join(tmp_dir, f"{environment_name}.yml")
            with open(yml_file, "w") as f:
                f.write(raw_spec)

            if self.environment_manager == "conda":
                return f"conda env create -f {yml_file} -n {environment_name}"

            elif self.environment_manager == "mamba":
                return f"mamba env create -f {yml_file} -n {environment_name}"

            elif self.environment_manager == "pixi":
                # Create a pixi project in the current working directory
                project_dir = os.path.join(os.getcwd(), environment_name)
                os.makedirs(project_dir, exist_ok=False)  # Raise exception if the project directory exists
                return (
                    f"pixi init --import {yml_file} {project_dir} && "
                    f"pixi install --manifest-path {project_dir}"
                )

        elif self.environment_manager == "uv":
            # Write the requirements.txt file
            req_file = os.path.join(tmp_dir, f"{environment_name}.txt")
            with open(req_file, "w") as f:
                f.write(raw_spec)
            # Create a uv project in the current working directory
            env_dir = os.path.join(os.getcwd(), environment_name)
            os.makedirs(env_dir, exist_ok=False)  # Raise exception if the project directory exists
            return (
                f"uv venv create {env_dir} && "
                f"uv pip sync -r {req_file} --venv {env_dir}"
            )

        raise RuntimeError(f"Unsupported environment manager: '{self.environment_manager}'")


@contextmanager
def not_on_worker() -> Generator[None, Any, None]:
    try:
        get_worker()
    except BaseException:
        yield

with not_on_worker():
    _env_config = EnvironmentConfig()
    environment_manager: str = _env_config.environment_manager
    environment_cmd: str = _env_config.env_export_cmd

source_domains: List[str] = [
    "conda.graylab.jhu.edu",
    "west.rosettacommons.org",
    "conda.rosettacommons.org",
]  # Conda channels and/or source domains (containing PyRosetta usernames/passwords) to be stripped from YML file strings.
