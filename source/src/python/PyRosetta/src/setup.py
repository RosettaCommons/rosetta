import os, sys, os.path, json, subprocess

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


def get_package_version():
    version_info = None

    package_version_file = os.path.dirname( os.path.abspath(__file__) ) + '/../version.json'
    if os.path.isfile(package_version_file):
        with open(package_version_file) as f: version_info = json.load(f)

    try:
        root_version_file = os.path.join(
                subprocess.check_output(["git", "rev-parse", "--show-toplevel"]).strip().decode(),
            "source/.version.json")
        if os.path.isfile(root_version_file):
            with open(root_version_file) as f: version_info = json.load(f)
    except subprocess.CalledProcessError:
        pass

    if version_info and version_info["version"]:
        return version_info['version']

    return 'v2017'


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
            version = get_package_version(),
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
