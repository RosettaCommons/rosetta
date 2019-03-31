from __future__ import print_function
import sys

if tuple(sys.version_info) < (3, 5):
    print("Unsupported python version for pyrosetta.distributed: %s" % sys.version_info)
    sys.exit(0)
else:
    print("sys.version_info: %s", sys.version_info)

import pyrosetta.rosetta
import shlex
import subprocess
import tempfile
import venv


if not hasattr(pyrosetta.rosetta, "cereal"):
    print("Unsupported non-serialization build for pyrosetta.distributed.")
    sys.exit(0)

def e(cmd):
    print(" ".join(map(shlex.quote, cmd)))
    subprocess.check_call(cmd)

test_suites = [
                "pyrosetta.tests.distributed.test_smoke",
                "pyrosetta.tests.distributed.test_dask",
                "pyrosetta.tests.distributed.test_concurrency",
                "pyrosetta.tests.distributed.test_gil",
                "pyrosetta.tests.bindings.core.test_pose",
                "pyrosetta.tests.numeric.test_alignment"
              ]

with tempfile.TemporaryDirectory(prefix="tmp_pyrosetta_env") as venv_dir:
    venv.create(venv_dir, clear=True, system_site_packages=True, with_pip=True)
    e(
        [
            venv_dir + "/bin/pip", "install",
            "blosc", "dask", "distributed", "numpy", "pandas", "scipy", "traitlets"
        ]
    )
    
    for test_suite in test_suites:
        e(
            [
                venv_dir + "/bin/python",
                "-m",
                "unittest",
                test_suite,
            ]
        )
