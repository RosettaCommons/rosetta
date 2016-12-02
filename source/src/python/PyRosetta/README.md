PyRosetta 4
===========

Building PyRosetta under Anaconda python
----------------------------------------
If you use Anaconda python (https://www.continuum.io) and would like to use
it to build PyRosetta, you will need to make the build system aware of the
location of the appropriate C headers and library. These values are specified
using the `python-include-dir` and `python-lib` command line flags when calling
`build.py`. For example, if Anaconda 3 is installed in the default location,
the command to run `build.py` would be:

`python build.py -jX --python-include-dir=$HOME/anaconda3/include/python3.5m --python-lib=$HOME/anaconda3/lib/libpython3.5m.dylib`

where X is the number of jobs to parallelize the build over.
