# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import logging
import os
import shutil
import sys
import warnings

from functools import lru_cache
from pyrosetta.utility import get_package_version
from typing import (
    Dict,
    Generic,
    List,
    Tuple,
    TypeVar,
)


__dask_version__: Tuple[int, int, int] = get_package_version("dask")
__dask_jobqueue_version__: Tuple[int, int, int] = get_package_version("dask-jobqueue")


G = TypeVar("G")


class EnvironmentConfig(Generic[G]):
    _ENV_VAR: str = "PYROSETTACLUSTER_ENVIRONMENT_MANAGER"
    _ENV_MANAGERS: Tuple[str, ...] = ("pixi", "uv", "mamba", "conda")
    _ENV_EXPORT_CMDS: Dict[str, str] = {
        "pixi": "pixi lock --check || pixi lock --no-install",
        "uv": "uv export --format requirements-txt --frozen",
        "mamba": f"mamba env export --prefix '{sys.prefix}'",
        "conda": f"conda env export --prefix '{sys.prefix}'",
    }

    def __init__(self) -> None:
        _env_var_manager = os.environ.get(EnvironmentConfig._ENV_VAR, None)
        if _env_var_manager:
            self.environment_manager = _env_var_manager
            logging.debug(
                "Configuring environment manager for PyRosettaCluster from operating system "
                + f"environment variable: {EnvironmentConfig._ENV_VAR}={self.environment_manager}"
            )
            if self.environment_manager not in EnvironmentConfig._ENV_MANAGERS:
                raise ValueError(
                    "The '{0}' environment variable must be in: '{1}'. Received: '{2}'.".format(
                        EnvironmentConfig._ENV_VAR,
                        EnvironmentConfig._ENV_MANAGERS,
                        self.environment_manager,
                    )
                )
        else:
            for _manager in EnvironmentConfig._ENV_MANAGERS:
                if shutil.which(_manager):
                    self.environment_manager = _manager
                    logging.debug(f"Configuring environment manager for PyRosettaCluster: '{_manager}'")
                    break
            else:
                self.environment_manager = "conda"
                warnings.warn(
                    f"Warning: could not configure an environment manager for PyRosettaCluster. "
                    + "Please ensure that either of 'pixi', 'uv', 'mamba', or 'conda' is installed. "
                    + "Using 'conda' as the default environment manager.",
                    UserWarning,
                    stacklevel=7,
                )

    @property
    def env_export_cmd(self) -> str:
        """
        Return the appropriate environment export command for the given environment manager.
        This method automatically adjusts for pixi and uv when a manifest path or project path
        is set via environment variables.

        Args:
            env_manager: A `str` object of either "pixi", "uv", "mamba", "conda".

        Returns:
            A shell command string for subprocess execution.
        """
        # Update pixi environment command if `$PIXI_PROJECT_MANIFEST` is set
        if self.environment_manager == "pixi":
            # https://pixi.sh/dev/reference/environment_variables/#environment-variables-set-by-pixi
            manifest_path = os.environ.get("PIXI_PROJECT_MANIFEST")
            if manifest_path:
                # Append `--manifest-path` flag to both commands in the OR clause
                logging.info(
                    "PyRosettaCluster detected the set 'PIXI_PROJECT_MANIFEST' environment variable, and is "
                    + f"setting the flag `--manifest-path '{manifest_path}'` in the `pixi lock` command."
                )
                return (
                    f"pixi lock --check --manifest-path '{manifest_path}' || "
                    f"pixi lock --no-install --manifest-path '{manifest_path}'"
                )

        # Update uv environment command if `$UV_PROJECT` is set
        elif self.environment_manager == "uv":
            # https://docs.astral.sh/uv/reference/environment/#uv_project
            project_dir = os.environ.get("UV_PROJECT")
            if project_dir:
                # Append `--project` flag
                logging.info(
                    "PyRosettaCluster detected the set 'UV_PROJECT' environment variable, and is "
                    + f"setting the flag `--project '{project_dir}'` in the `uv export` command."
                )
                return (
                    f"uv export --format requirements-txt --frozen --project '{project_dir}'"
                )

        # Use default environment export command
        return self._ENV_EXPORT_CMDS[self.environment_manager]


@lru_cache(maxsize=1)
def get_environment_config() -> EnvironmentConfig:
    """Return an instance of the `EnvironmentConfig` class on the host process."""
    return EnvironmentConfig()


def get_environment_manager() -> str:
    """Get the configured environment manager."""
    return get_environment_config().environment_manager


def get_environment_cmd() -> str:
    """Get the configured environment export command."""
    return get_environment_config().env_export_cmd


def get_environment_var() -> str:
    """Get the PyRosettaCluster operating system environment variable name."""
    return EnvironmentConfig._ENV_VAR


source_domains: List[str] = [
    "conda.graylab.jhu.edu",
    "west.rosettacommons.org",
    "conda.rosettacommons.org",
]  # Conda channels and/or source domains (potentially containing PyRosetta usernames/passwords) to be sanitized from environment file strings.
