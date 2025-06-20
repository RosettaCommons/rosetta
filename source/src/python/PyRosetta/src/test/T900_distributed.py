from __future__ import print_function
import shutil
import sys


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


# check if required packages was installed
try:
    import attrs
    import cloudpickle
    import dask
    import pandas
except ImportError as e:
    print(f"{e}\nSome packages required for pyrosetta.distributed is missing, skipping the tests...")
    sys.exit(0)


if shutil.which("conda"):
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


def e(cmd):
    print("Executing:\n{0}".format(cmd))
    status, output = subprocess.getstatusoutput(cmd)
    print("Output:\n{0}".format(output))
    if status != 0:
        print(
            "Encountered error(s) with exit code {0} while running: {1}\nTerminating...".format(
                status, cmd
            )
        )
        sys.exit(1)


test_suites = [
    "pyrosetta.tests.bindings.core.test_pose",
    "pyrosetta.tests.distributed.cluster.test_logging",
    "pyrosetta.tests.distributed.cluster.test_reproducibility",
    "pyrosetta.tests.distributed.cluster.test_smoke",
    "pyrosetta.tests.distributed.test_concurrency",
    "pyrosetta.tests.distributed.test_dask",
    "pyrosetta.tests.distributed.test_gil",
    "pyrosetta.tests.distributed.test_smoke",
    "pyrosetta.tests.distributed.test_viewer",
    "pyrosetta.tests.numeric.test_alignment",
]


for test_suite in test_suites:
    e("{python} -m unittest {test_suite}".format(python=sys.executable, test_suite=test_suite))
    #import unittest
    #unittest.main(test_suite)
