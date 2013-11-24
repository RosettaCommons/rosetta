import warnings

from setuptools import setup, find_packages
import subprocess
import datetime

import os
from os import path

# Detect version number
MAJOR = 3
MINOR = 4
ISRELEASED = False

def resolve_build_version():
    from rosetta.version import commit
    version = commit

    if subprocess.call("git rev-parse --git-dir", shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE) == 0:
        import dateutil
        import dateutil.parser
        # If rev-parse --git-dir does not return an error we are in a git directory,
        # use rev-parse to get short version number and resolve the git commit date
        version = subprocess.check_output("git rev-parse --short %s" % version, shell=True).strip()
        head_version = subprocess.check_output("git rev-parse --short HEAD", shell=True).strip()

        if version != head_version:
            warnings.warn("Rosetta build version: %s does not match current git head: %s" % (version, head_version))

        commit_date = subprocess.check_output("git show -s --format='%%ci' %s" % version, shell=True).strip()
        commit_time = dateutil.parser.parse(commit_date).astimezone(dateutil.tz.tzutc())
        version_time = commit_time.strftime("%Y%m%d%H%M%S")
        full_version = "%s.%s.%s.git-%s" % (MAJOR, MINOR, version_time, version)
    else:
        warnings.warn("Unable to resolve version commit date, build in repository working directory.")
        full_version = "%s.%s.git-%s" % (MAJOR, MINOR, version)

    return full_version


if not ISRELEASED:
    full_version = resolve_build_version()
else:
    full_version = "%s.%s" % (MAJOR, MINOR)

# Extract list if all database files, can not write a simple glob
# that matches all db files but not db sub-directories.
database_files = []
for dirpath, _, files in os.walk("database"):
    dirrel = dirpath[len("database/"):]
    database_files.extend(path.join(dirrel, f) for f in files)

setup(
        name="pyrosetta",
        description="PyRosetta distributable package.",
        version = full_version,
        author_email="fordas@uw.edu",
        packages=find_packages() + ["database"],
        package_data = {
            "" : ["*.so", "libboost_python*"],
            "database" : database_files
            },
        zip_safe = False
        )
