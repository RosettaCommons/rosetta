#!/usr/bin/env python2.7

#Build Rosetta using CMake + Ninja.

import argparse
import subprocess
import sys
import os
import time

start_time = time.time()

class NinjaBuildError(Exception):
    pass

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def clean_build_path(build_path):
    for root, dirs, files in os.walk(build_path, topdown=False):
        for name in files:
            if name != "CMakeLists.txt" and name != "buildme.sh":
                os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))

def exe_exists(exe):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        exe_file = os.path.join(path, exe)
        if is_exe(exe_file):
            return True
    return False

#Input argument parsing with argparse
help_message1 = "Build Rosetta using CMake + Ninja. See source/cmake/README for more info."
help_message2 = (
    "examples:\n" +
    "ninja.build.py release        Standard Release build\n" +
    "ninja.build.py r              Same as above\n" +
    "ninja.build.py r -my          Release build with *.src.settings.my\n" +
    "ninja.build.py debug -remake  Debug build, rerun cmake + make_project.py\n" )
#Useful abbreviations
#Unit and debug are now combined
abbrev = {'r':'release', 'd':'debug', 'u':'debug'}

help_message2 += "\nbuild name abbreviations:\n"
for key, item in abbrev.items():
    help_message2 += '%s\t== %s\n' % (key, item)


parser = MyParser(description=help_message1, epilog=help_message2, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('build', default='release',
        help="Build name to be compiled in the cmake folder.")
parser.add_argument('-t', '-target', type=str, default=None,
        help="The build target to tell ninja to compile. (e.g. 'bin', 'unit', 'apps', 'pilot_apps', 'relax') Default is to compile everything.")
parser.add_argument('-remake', action='store_true', help="Rerun make_project.py and cmake. " +
        "Be sure to call it when building first time, after editing src.settings, or first time switching to -my.")
parser.add_argument('-my', action='store_true', help="Use instead *.src.settings.my in make_project.py.")
parser.add_argument('-clean', action='store_true', help="Remove all old files in build folder.")
parser.add_argument('-clean_exit', action='store_true', help="Remove all old files in build folder and exit. Overides all building options.")
parser.add_argument('-v', action='store_true', help="Where possible, run the verbose version of the commands.")
parser.add_argument('-j', type=int, default=None, help="The -j value to pass to ninja. NOTE: Ninja is good at autodiscovery for this value. Don't set it unless you know that's what you need.")
parser.add_argument('-k', action='store_true', help="Tell ninja to keep going as much as possible when it encounters an error, instead of stopping immediately. Good for lunch-break compiles.")
args = parser.parse_args()

#Convert abbreviations
if args.build in abbrev:
    args.build = abbrev[args.build]

#Check if required executable exists
required_executables = ('cmake', 'ninja')
for executable in required_executables:
    if not exe_exists(executable):
        raise NinjaBuildError("Cannot find executable %s! Ensure it is properly installed and added to your PATH." % executable)

#Check if required files exists
source_path = os.path.split(os.path.abspath(__file__))[0]
cmake_path = os.path.join(source_path, "cmake")
build_path = os.path.join(cmake_path, "build_%s" % args.build)

required_files = ( os.path.join(cmake_path, 'make_project.py'), os.path.join(build_path, 'CMakeLists.txt') )
for filename in required_files:
    if not os.path.isfile(filename):
        raise NinjaBuildError("Cannot find required file %s! Check if your Rosetta folder is complete." % filename)

#Clean the build path folder if requested
if args.clean_exit or args.clean:
    clean_build_path(build_path)
    if args.clean_exit:
        sys.exit(0)

#Run make_project.py if requested
if args.remake:
    os.chdir(cmake_path)
    if args.my:
        subprocess.check_call(["./make_project.py", "my"])
    else:
        subprocess.check_call(["./make_project.py", "all"])

#Run cmake if requested
os.chdir(build_path)
if args.remake:
    subprocess.check_call(["cmake", "-G", "Ninja"])
elif not os.path.isfile("build.ninja"):
    raise NinjaBuildError("File build.ninja not found. Use option '-remake' the first time building.")

#Run ninja to build Rosetta
ninja_command = ["ninja"]

if args.j is not None:
    ninja_command.extend( ["-j",str(args.j)] )

if args.t is not None:
    ninja_command.append( str(args.t) )

if args.v:
    ninja_command.append( "-v" )
if args.k:
    ninja_command.extend( ["-k","5000"] )

subprocess.check_call(ninja_command)

print "Job completed. Total time = %f s" % ( time.time() - start_time )

