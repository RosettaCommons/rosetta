import os, sys

# Assert or bootstrap minimum setuptools version required for find_packages
import ez_setup
ez_setup.use_setuptools("3.3")

from setuptools import setup, find_packages, Distribution

def get_file_list(directory):
    """Get list of all files in directory, relative to directory name."""
    result_files = []
    for dirpath, _, files in os.walk(directory):
        dirrel = dirpath[len(directory) + 1:]
        result_files.extend(os.path.join(dirrel, f) for f in files)

    return result_files

class PyRosettaDistribution( Distribution ):
    """Override default distribution class to ensure packages are marked as containing non-pure modules."""
    def is_pure(self):
        return False

    def has_ext_modules(self):
        return True

def setup_package():
    # chdir to source path and add sources to sys.path to allow invocation of setup script from other directories.
    # Code structure taken from numpy/setup.py
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    try:
        setup(
            name = "pyrosetta",
            description="PyRosetta package",
            version = "4.0",
            packages = find_packages(include=["pyrosetta*", "rosetta*"]),
            package_data = {
                "pyrosetta" :
                    ["*.so"] +
                    [os.path.join("database", f) for f in get_file_list("pyrosetta/database")] +
                    [os.path.join("lib", f) for f in get_file_list("pyrosetta/lib")],
            },
            distclass = PyRosettaDistribution,
            zip_safe = False,
        )
    finally:
        os.chdir(old_path)
        del sys.path[0]

if __name__ == "__main__":
    setup_package()
