
Rosetta/PyRosetta reference Docker recipes
==========================================

Collection of various recipes to serve as example for setting up Rosetta/PyRosetta build environment.
Note that these recipes _does not actually clone/build Rosetta or PyRosetta_ and _only serve as examples_ to
how-to setup build environments. When using this images you will need to mount already cloned Rosetta repository build commands:
```
# to build Rosetta
cd rosetta/source && ./scons.py -j8 mode=release bin

# to build PyRosetta
cd rosetta/source/src/python/PyRosetta && python3 build.py -j8
```

* `alpine-3.9` absolutely minimal Rosetta build environment
* `ubuntu-xx.xx` reference images for building Rosetta and PyRosetta
* `ml` PyRosetta reference image to build PyRosetta with `libtorch` and `TensorFlow`  support
