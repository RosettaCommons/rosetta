from __future__ import print_function
import os
import shutil
import sys
import tempfile
import time


if tuple(sys.version_info) < (3, 6):
    print("Unsupported python version for pyrosetta.distributed: %s" % str( sys.version_info ) )
    sys.exit(0)
else:
    print("sys.version_info: {0}".format(sys.version_info))

import pyrosetta.rosetta
import subprocess

if not hasattr(pyrosetta.rosetta, "cereal"):
    print("Unsupported non-serialization build for pyrosetta.distributed.")
    sys.exit(0)

try:
    import attrs
    import cloudpickle
    import dask
    import pandas
except ImportError as e:
    print(f"{e}\nSome packages required for pyrosetta.distributed is missing, skipping the tests...")
    sys.exit(0)

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
                print("Current pixi environment:", export, sep="\n")
            else:
                print("The exported pixi environment file does not exist: '{0}'.".format(tmp_file))
    except subprocess.CalledProcessError as ex:
        print("Printing pixi environment failed with return code: {0}.".format(ex.returncode))
elif shutil.which("conda"):
    try:
        export = subprocess.check_output(
            "conda env export --prefix $(conda env list | grep '*' | awk '{print $NF}')",
            shell=True,
            text=True,
        )
        print("Current conda environment:", export, sep="\n")
    except subprocess.CalledProcessError as ex:
        print("Printing conda environment failed with return code: {0}.".format(ex.returncode))
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


tests = [
    "pyrosetta.tests.billiard.stress_test"
    "pyrosetta.tests.billiard.stress_test_2"
]
for test in tests:
    t0 = time.perf_counter()
    run_subprocess("{python} -m unittest {test}".format(python=sys.executable, test=test))
    t1 = time.perf_counter()
    dt = t1 - t0
    print("Finished running test in {0} seconds: {1}\n".format(round(dt, 6), test))
