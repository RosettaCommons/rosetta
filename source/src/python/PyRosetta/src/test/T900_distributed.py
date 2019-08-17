from __future__ import print_function
import sys

if tuple(sys.version_info) < (3, 5):
    print("Unsupported python version for pyrosetta.distributed: %s" % sys.version_info)
    sys.exit(0)
else:
    print("sys.version_info: %s", sys.version_info)

import pyrosetta.rosetta
import subprocess
import tempfile
import venv

if not hasattr(pyrosetta.rosetta, "cereal"):
    print("Unsupported non-serialization build for pyrosetta.distributed.")
    sys.exit(0)

def e(cmd):
    print("Executing:\n{0}".format(cmd))
    status, output = subprocess.getstatusoutput(cmd)
    print("Output:\n{0}".format(output))
    if status != 0:
        print("Encountered error(s) with exit code {0} while running: {1}\nTerminating...".format(status, cmd))
        sys.exit(1)


test_suites = [
                "pyrosetta.tests.bindings.core.test_pose",
                "pyrosetta.tests.distributed.test_concurrency",
                "pyrosetta.tests.distributed.test_dask",
                "pyrosetta.tests.distributed.test_gil",
                "pyrosetta.tests.distributed.test_smoke",
                "pyrosetta.tests.distributed.test_viewer",
                "pyrosetta.tests.numeric.test_alignment"
              ]

with tempfile.TemporaryDirectory(prefix="tmp_pyrosetta_env") as venv_dir:
    
    venv.create(venv_dir, clear=True, system_site_packages=False, with_pip=True)
    
    packages = "blosc dask distributed jupyter numpy pandas py3Dmol scipy traitlets"
    e(
        "source {0}/bin/activate && {0}/bin/pip install {1}".format(venv_dir, packages)
    )
    
    for test_suite in test_suites:
        e(
            "source {0}/bin/activate && {0}/bin/python -m unittest {1}".format(venv_dir, test_suite)
        )
