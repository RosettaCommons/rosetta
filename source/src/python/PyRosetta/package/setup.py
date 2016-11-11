import os, sys
from setuptools import setup, find_packages, Distribution, Extension


# Check if current working dir is where this file is located...
if os.path.dirname(os.path.realpath(__file__)) != os.getcwd():
    print('This script must be run from within setup directory! Exiting...')
    sys.exit(1)


def get_file_list(directory):
    """Get list of all files in directory, relative to directory name."""
    result_files = []
    for dirpath, _, files in os.walk(directory):
        dirrel = dirpath[len(directory) + 1:]

        dir_ = os.path.join( os.path.basename(directory), os.path.join(dirrel) )
        result_files.append( (dir_, [os.path.join(dir_, f) for f in files]) )

    return result_files


class PyRosettaDistribution( Distribution ):
    """Override default distribution class to ensure packages are marked as containing non-pure modules."""
    def is_pure(self): return False


setup(
    name = "pyrosetta",
    description="PyRosetta package",
    version = "4.0",
    packages = find_packages(),
    #scripts = ['rosetta.pyd'],
    #py_modules=['rosetta'],

    data_files = [('.', ['rosetta.so'])] + get_file_list("database"),
    # data_files=[#('.', ['rosetta.so'])
    #             #('.', [os.path.join("database", f) for f in get_file_list("database")]), ],
    #data_files=[('database/gpu', ['database/gpu/DARC_PSO.cl', 'database/gpu/darc.cl',])],

    # ext_modules=[
    #     Extension(name="rosetta",
    #               sources=[],
    #               #include_dirs=['.'],
    #               #define_macros=DEFINE_MACROS,
    #               #extra_objects=['readline/libreadline.a', 'readline/libhistory.a'],
    #               #libraries=['ncurses']
    #     ),
    # ],
    package_data = {
         #"" :    ["*.so"] + [f for f in get_file_list("database")]
         #"" : [os.path.join("database", f) for f in get_file_list("database")],
    },
    distclass = PyRosettaDistribution,
    zip_safe = False,
    )
