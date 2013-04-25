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
print "unit: build unit test version."
print "###################################################################"

#Command line options
is_debug = False
is_remake = False
is_unit = False
if "remake" in sys.argv :
    is_remake = True
if "debug" in sys.argv:
    is_debug = True
if "unit" in sys.argv:
    is_unit = True

#Check if required files exists
assert exists("./cmake/make_project.py")
assert exists("./cmake/build_release/CMakeLists.txt")
assert exists("./cmake/build_debug/CMakeLists.txt")
assert exists("./cmake/build_unit/CMakeLists.txt")

os.chdir("./cmake")
if is_remake :
    subprocess.check_call("./make_project.py all", shell=True)

if is_unit:
    os.chdir("./build_unit")
elif is_debug :
    os.chdir("./build_debug")
else :
    os.chdir("./build_release")

if is_remake :
    subprocess.check_call("cmake -G Ninja", shell=True)

if not os.path.exists("build.ninja"):
    print "------------------------------------------------------------------------------------------------------"
    print "----- File build.ninja not found: Call with the 'remake' option the first time running a build. ------"
    print "------------------------------------------------------------------------------------------------------"
    sys.exit(1)

subprocess.check_call("ninja")

print "DONE!...Total time = %f s" % ( time.time() - start_time )
