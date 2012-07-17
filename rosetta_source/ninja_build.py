#!/usr/bin/env python

#Build Rosetta using CMake + Ninja.

import subprocess
import sys
from os.path import exists
import os
import time

start_time = time.time()

print "###################################################################"
print "Build Rosetta using CMake + Ninja."
print "See rosetta_source/cmake/README_ninja for setup instruction."
print "Default - run ninja only to build release version without reloading *.src.settings"
print "Options"
print "remake: Full build process by running make_project.py and cmake too (slower)."
print "debug:  build debug version instead."
print "###################################################################"

#Command line options
is_debug = False
is_remake = False
if "remake" in sys.argv :
    is_remake = True
if "debug" in sys.argv:
    is_debug = True

#Check if required files exists
assert exists("./cmake/make_project.py")
assert exists("./cmake/build_release/CMakeLists.txt")
assert exists("./cmake/build_debug/CMakeLists.txt")

os.chdir("./cmake")
if is_remake :
    subprocess.check_call("./make_project.py all", shell=True)

if is_debug :
    os.chdir("./build_debug")
else :
    os.chdir("./build_release")

if is_remake :
    subprocess.check_call("cmake -G Ninja", shell=True)

subprocess.check_call("ninja")

print "DONE!...Total time = %f s" % ( time.time() - start_time )
