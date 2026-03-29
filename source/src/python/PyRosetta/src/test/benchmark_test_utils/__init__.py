# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import os
import shutil
import subprocess
import sys
import tempfile
import time

from pyrosetta.utility import has_cereal
from typing import NoReturn, Optional


def has_pyrosetta_distributed_package_requirements() -> bool:
    """Test if `pyrosetta.distributed` framework packages are installed in the virtual environment."""
    try:
        import attrs
        import billiard
        import blosc
        import cloudpickle
        # import cryptography
        import dask
        import dask_jobqueue
        import distributed
        import git
        import jupyter
        import numpy
        import pandas
        # import py3Dmol
        import scipy
        import traitlets
        # import lzma as xz
        return True
    except ImportError as ex:
        print(f"{type(ex).__name__}: {ex}")
        return False


def has_python_version(major: int, minor: int) -> bool:
    """Test if the Python version is greater than or equal to the provided `major` and `minor` values."""
    return tuple(sys.version_info) >= (major, minor)


def exit_if_missing_pyrosetta_distributed_requirements(returncode: int = 0) -> Optional[NoReturn]:
    """Exit the Python process if the `pyrosetta.distributed` framework requirements are missing from the virtual environment."""

    if not has_python_version(3, 6):
        print("Unsupported Python version for `pyrosetta.distributed` framework: %s" % str(sys.version_info))
        sys.exit(returncode)
    if not has_cereal():
        print("Unsupported non-serialization build for `pyrosetta.distributed` framework test.")
        sys.exit(returncode)
    if not has_pyrosetta_distributed_package_requirements():
        print(f"Packages required for `pyrosetta.distributed` framework tests are missing. Skipping the tests...")
        sys.exit(returncode)


def print_environment_export() -> None:
    if shutil.which("pixi"):
        try:
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_file = os.path.join(tmp_dir, "environment.yml")
                subprocess.run(
                    f"pixi workspace export conda-environment '{tmp_file}'",
                    shell=True,
                )
                if os.path.isfile(tmp_file):
                    with open(tmp_file, "r") as f:
                        export = f.read()
                    print("Current Pixi environment:", export, sep="\n")
                else:
                    print("The exported Pixi environment file does not exist: '{0}'.".format(tmp_file))
        except subprocess.CalledProcessError as ex:
            print("Printing Pixi environment failed with return code: {0}.".format(ex.returncode))
    elif shutil.which("conda"):
        try:
            export = subprocess.check_output(
                "conda env export --prefix $(conda env list | grep '*' | awk '{print $NF}')",
                shell=True,
                text=True,
            )
            print("Current Conda environment:", export, sep="\n")
        except subprocess.CalledProcessError as ex:
            print("Printing Conda environment failed with return code: {0}.".format(ex.returncode))
    else:
        try:
            freeze = subprocess.check_output(
                "{0} -m pip freeze".format(sys.executable),
                shell=True,
                text=True,
            )
            print("Current pip environment:", freeze, sep="\n")
        except subprocess.CalledProcessError as ex:
            print("Printing pip environment failed with return code: {0}.".format(ex.returncode))


def run_subprocess(cmd: str) -> int:
    """Run a command in a subprocess."""

    print("Executing:\n{0}".format(cmd))
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    print("Output:")
    for line in process.stdout:
        print(line, end="")

    process.wait()
    status = process.returncode
    if status != 0:
        print(f"Encountered error(s) with exit code {status} while running: {cmd}\nTerminating...")
        sys.exit(1)


def run_test_cases(*test_cases: str) -> None:
    """Run the input test cases using the `unittest` module."""

    exit_if_missing_pyrosetta_distributed_requirements()
    print_environment_export()

    for test_case in test_cases:
        t0 = time.perf_counter()
        run_subprocess(f"{sys.executable} -m unittest {test_case}")
        t1 = time.perf_counter()
        dt = t1 - t0
        print(f"Finished running test in {dt:.6f} seconds: {test_case}\n")


def run_distributed_cluster_test_cases(*test_cases: str) -> None:
    """Run the input test cases (each prepended with "pyrosetta.tests.distributed.cluster.") using the `unittest` module."""

    prefix = "pyrosetta.tests.distributed.cluster."
    prepend = lambda test_case: f"{prefix}{test_case}"
    run_test_cases(*map(prepend, test_cases))
