PyRosetta 4
===========


**NOTE: Use for PyRosetta for commercial purposes requires purchase of a separate PyRosetta license which is different from Rosetta commercial license.** (This includes fee-for-service work by academic users.) Please see [UW Comotion](https://els2.comotion.uw.edu/product/pyrosetta) or email license@uw.edu for more information.

Building from source
--------------------

Requirements:
* CLang-3.4+ or GCC-4.8+
* CMake
* [Ninja](https://ninja-build.org/)
* Python-3.8+

Building PyRosetta:
```
cd rosetta/source/src/python/PyRosetta
python3 build.py -j8
```

Useful build script options:
* `--create-package` create Python setup.py package
* `--create-wheel` create Python wheel package
* `--print-build-root` print path to where PyRosetta binaries will be located with given build options


Developing PyRosetta
---------------------------
When developing/debugging PyRosetta is usually too time consuming to create and install it as a package during testing. To mitigate this PyRosetta build system is setup in way to allow PyRosetta import to be done at build location. Also, all Python files at build location are symlinks to original source files so it is possible to edit and test Python scripts without explicitly rebuilding the PyRosetta.

To use this feature you will first need to determent path of binaries with currently build options using `--print-build-root` build script options. Note that report path will be vary for different build options, compiler settings etc. Final build path could be obtained by adding `/build` suffix to reported location. Starting python interpreter at final build path will allow you to fully import and PyRosetta as if it was installed.
