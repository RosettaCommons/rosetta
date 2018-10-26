from __future__ import print_function
import sys
import os

if tuple(sys.version_info) < (3, 5):
    print("Unsupported python version for pyrosetta.distributed: %s" % sys.version_info)
    sys.exit(0)
else:
    print("sys.version_info: %s", sys.version_info)

import tempfile
import venv
import subprocess
import shlex

import pyrosetta.rosetta

if not hasattr(pyrosetta.rosetta, "cereal"):
    print("Unsupported non-serialization build for pyrosetta.distributed.")
    sys.exit(0)


def e(cmd):
    print(" ".join(map(shlex.quote, cmd)))
    subprocess.check_call(cmd)


with tempfile.TemporaryDirectory(prefix="tmp_pyrosetta_env") as venv_dir:
    venv.create(venv_dir, clear=True, system_site_packages=True, with_pip=True)
    e([venv_dir + "/bin/pip", "install", "numpy", "pandas", "traitlets"])

    e(
        [
            venv_dir + "/bin/python",
            "-m",
            "unittest",
            "pyrosetta.tests.distributed.test_smoke",
        ]
    )
