# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import functools
import os
import select
import shutil
import signal
import subprocess
import sys
import tempfile
import time

from io import TextIOWrapper
from pyrosetta.utility import has_cereal
from typing import (
    Any,
    Callable,
    List,
    NoReturn,
    Optional,
    TypeVar,
    cast,
)

F = TypeVar("F", bound=Callable[..., int])

TIMEOUT: int = 1800 # seconds
STREAMING: bool = False


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
    """Print the active virtual environment export."""

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


def handle_status(func: F) -> F:
    """Decorator for test timing and status reporting."""

    @functools.wraps(func)
    def wrapper(test_case: str, *args: Any, **kwargs: Any) -> int:
        """Wrapper function for `handle_status` decorator."""

        t0 = time.perf_counter()
        status = func(test_case, *args, **kwargs)
        t1 = time.perf_counter()
        dt = t1 - t0

        if status == 0:
            print(f"Finished running test in {dt:.6f} seconds: {test_case}\n", flush=True)
        else:
            print(
                f"Encountered error(s) with exit code {status} after {dt:.6f} seconds while running: {test_case}",
                "Terminating...",
                sep="\n",
                flush=True,
            )
            sys.exit(1)

        return status

    return cast(F, wrapper)


def terminate_process(pid: int) -> None:
    """Terminate a process group."""

    try:
        os.killpg(os.getpgid(pid), signal.SIGTERM)
        time.sleep(1)
        os.killpg(os.getpgid(pid), signal.SIGKILL)
    except Exception:
        pass


def drain_buffer(pipe: TextIOWrapper) -> None:
    """Drain buffer from a pipe."""
    try:
        while True:
            rlist, _, _ = select.select([pipe], [], [], 0)
            if not rlist:
                break
            line = pipe.readline()
            if not line:
                break
            print(line, end="", flush=True)
    except Exception:
        pass


@handle_status
def run_unittest_getstatusoutput(test_case: str) -> int:
    """
    Run a test case using the `unittest` module with `subprocess.getstatusoutput`.

    Args:
        `test_case`: `str`
            The test module to run.

    Returns:
        An `int` object representing the return code.
    """

    args: List[str] = [sys.executable, "-m", "unittest", test_case]
    cmd: str = " ".join(args)
    print("Executing:", cmd, sep="\n", flush=True)

    status, output = subprocess.getstatusoutput(cmd)

    print("Output:", flush=True)
    print(output, end="", flush=True)

    return status


@handle_status
def run_unittest(test_case: str, timeout: int) -> int:
    """
    Run a test case using the `unittest` module in a subprocess.

    Args:
        `test_case`: `str`
            The test module to run.

        `timeout`: `int`
            The timeout in seconds.

    Returns:
        An `int` object representing the return code.
    """

    args: List[str] = [sys.executable, "-m", "unittest", test_case]
    print("Executing:", " ".join(args), sep="\n", flush=True)

    process = subprocess.Popen(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        start_new_session=True,
        close_fds=True,
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
    )

    try:
        try:
            stdout, _ = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            print(
                f"Subprocess timeout after {timeout} seconds while running: {test_case}",
                "Terminating...",
                sep="\n",
                flush=True,
            )
            terminate_process(process.pid)
            stdout, _ = process.communicate()
            print(stdout, end="", flush=True)
            return 1

        print("Output:", flush=True)
        print(stdout, end="", flush=True)

    finally: # Clean up in case any subprocesses are still alive
        if process.poll() is None:
            terminate_process(process.pid)

    return process.returncode


@handle_status
def run_unittest_streaming(test_case: str, timeout: int) -> int:
    """
    Run a test case using the `unittest` module in a subprocess with streaming standard output.

    Args:
        `test_case`: `str`
            The test module to run.

        `timeout`: `int`
            The timeout in seconds.

    Returns:
        An `int` object representing the return code.
    """

    args: List[str] = [sys.executable, "-m", "unittest", test_case]
    print("Executing:", " ".join(args), sep="\n", flush=True)

    process = subprocess.Popen(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        start_new_session=True,
        bufsize=1,
        close_fds=True,
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
    )
    if process.stdout is None:
        return process.wait()

    start_time = time.perf_counter()
    try:
        print("Output:", flush=True)
        while True:
            # Wait for output readiness without blocking
            rlist, _, _ = select.select([process.stdout], [], [], 0.1)
            if rlist:
                line = process.stdout.readline()
                if line:
                    print(line, end="", flush=True)

            # Check timeout
            if time.perf_counter() - start_time > timeout:
                print(
                    f"Subprocess timeout after {timeout} seconds while running: {test_case}",
                    "Terminating...",
                    sep="\n",
                    flush=True,
                )
                terminate_process(process.pid)
                time.sleep(0.1)
                drain_buffer(process.stdout) # Drain remaining buffer
                return 1

            # Check if process has exited
            if process.poll() is not None:
                break

        drain_buffer(process.stdout) # Drain remaining buffer after exit

    finally: # Clean up in case any subprocesses are still alive
        if process.poll() is None:
            terminate_process(process.pid)

    return process.returncode


def run_test_cases(*test_cases: str, streaming: bool = STREAMING, timeout: int = TIMEOUT) -> None:
    """
    Run the input test cases using the `unittest` module.

    Args:
        `*test_cases`: `str`
            The test modules to run.

        `streaming`: `bool`
            Whether or not to run the test modules with streaming standard output.

        `timeout`: `int`
            The timeout in seconds for each test module.

    Returns:
        `None`
    """

    exit_if_missing_pyrosetta_distributed_requirements()
    print_environment_export()
    if streaming:
        for test_case in test_cases:
            run_unittest_streaming(test_case, timeout)
    else:
        for test_case in test_cases:
            run_unittest_getstatusoutput(test_case)
    sys.stdout.flush()


def run_distributed_cluster_test_cases(*test_cases: str, streaming: bool = STREAMING, timeout: int = TIMEOUT) -> None:
    """
    Run the input test cases (each prepended with "pyrosetta.tests.distributed.cluster.") using the `unittest` module.

    Args:
        `*test_cases`: `str`
            The test modules to run (each prepended with "pyrosetta.tests.distributed.cluster.")

        `streaming`: `bool`
            Whether or not to run the test modules with streaming standard output.

        `timeout`: `int`
            The timeout in seconds for each test module.

    Returns:
        `None`
    """

    prefix = "pyrosetta.tests.distributed.cluster."
    run_test_cases(
        *(f"{prefix}{test_case}" for test_case in test_cases),
        streaming=streaming,
        timeout=timeout,
    )
