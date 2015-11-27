#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true: :collapseFolds=10:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   BuildBuindings.py
## @brief  Build Python buidings for rosetta
## @author Sergey Lyskov

import logging
try:
    import coloredlogs
    log_config = coloredlogs.install(
        level=logging.INFO,
        fmt="[%(process)d] %(asctime)s %(filename)s:%(lineno)s %(levelname)s %(message)s"
    )
except ImportError:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(process)d] %(asctime)s %(filename)s:%(lineno)s %(levelname)s %(message)s"
    )

import os, re, sys, time, commands, shutil, platform, os.path, gc, json, glob
import pprint
from os import path
import subprocess #, errno

import multiprocessing
import signal

import cStringIO as StringIO
# See https://github.com/jreese/multiprocessing-keyboardinterrupt/
def init_worker():
    signal.signal( signal.SIGINT, signal.SIG_IGN )

# Expand path to include tools subdir
script_root_dir = path.dirname(path.realpath(__file__))
logging.info("Script root directory: %s", script_root_dir)
sys.path.append( path.dirname( script_root_dir  ))

# Create global 'current_platform' that will hold info of current system
if sys.platform.startswith("linux"):
    current_platform = "linux" # can be linux1, linux2, etc
    current_lib_suffix = "so"
elif sys.platform == "darwin" :
    current_platform = "macos"
    current_lib_suffix = "dylib"
elif sys.platform == "cygwin" :
    raise ValueError("Unsupportred platform: %s", sys.platform)
    current_platform = "cygwin"
    current_lib_suffix = 'dll'
elif sys.platform == "win32" :
    raise ValueError("Unsupportred platform: %s", sys.platform)
    current_platform = "windows"
    current_lib_suffix = 'dll'
else:
    raise ValueError("Unrecognized platform: %s", sys.platform)

current_platform_bits = platform.architecture()[0][:2]

import exclude

import tools.DoxygenExtractorPyPP
import tools.CppParser

from argparse import ArgumentParser

from collections import namedtuple
JobTuple = namedtuple("JobTuple",["pid", "tag"])

candidate_target_modules = ['utility', 'numeric', 'basic', 'core', 'protocols']
candidate_target_library = {
    'numeric' : "core.5",
    'basic' : "core.5",
    'core' : "core.5",
    'utility' : "core.5",
    'protocols' : "protocols.7"
}

def main(args):
    ''' Script to build Python buidings.
    '''
    parser = ArgumentParser(usage="Generate pyrosetta distribution.")

    parser.add_argument('-I',
      default=[],
      action="append",
      help="Additional include paths used in bindings build.",
    )

    parser.add_argument('-L',
      default=[],
      action="append",
      help="Additional libraries paths used in bindings link.",
    )

    parser.add_argument("-1", "--one",
      default=[], action="append",
      help="Build just one namespace instead of whole project, can be specified multiple times.",
    )

    parser.add_argument("-t", "--target",
        type=str, default=None, choices=candidate_target_modules,
        help="Target build namespace, dependent namespaces built. (Choice")

    parser.add_argument("--build_rosetta",
      default=True,
      action="store_true",
      help="Build rosetta libraries before generating bindings."
      )

    parser.add_argument("-d",
      action="store_false", dest="build_rosetta",
      help="Distable build rosetta libraries before generating bindings.",
    )

    parser.add_argument("-u", "--update",
      default=False,
      action="store_true",
      help="Debug only. Try to check time stamp of files before building them.",
      )

    parser.add_argument("--debug",
      action="store_true", default=False,
      help="Perform a Debug build when possible.",
      )

    parser.add_argument("--debug_bindings",
      action="store_true", default=False,
      help="Build bindings with -DDEBUG and debug symbols.",
      )

    parser.add_argument("--numpy_support",
      action="store_true", default=True,
      help="Enable numpy type conversion support.",
      )

    parser.add_argument("--no-numpy_support",
      action="store_false", default=False, dest="numpy_support",
      help="Disable numpy type conversion support.",
      )

    parser.add_argument("--package_path",
        action="store", default="pyrosetta",
        help="Target package directory for bindings. (default: %(default)s)",)

    parser.add_argument("--continue_on_error",
      default=False,
      action="store_true", dest="continue_on_error",
      help="Debug only. Continue building after compilation errors."
      )

    parser.add_argument("--compiler",
      default='gcc',
      action="store",
      help="Default compiler that will be used to build PyRosetta. (default: %(default)s)",
      )

    parser.add_argument(
            "--boost_path",
            action="store",
            required=False,
            default=None,
            help="Path to boost install prefix.")

    parser.add_argument("--boost_lib",
      default='boost_python',
      action="store",
      help="Name of boost dynamic library. (default: %(default)s)",
      )

    parser.add_argument(
            "--package_boost",
            action="store_true", default=False,
            help="Package external boost_python library specified under --boost_path."
    )

    parser.add_argument(
      "--python_lib",
      default="python2.7",
      action="store",
      help="Target python version. (default: %(default)s)",
    )

    parser.add_argument(
        "--python_path",
        action="store",
        required=False,
        default=None,
        help="Python install prefix.")

    parser.add_argument("--max-function-size", default=1024*128, type=int,
            help="Maximum size of function in binding files in bytes. (default: %(default)s)"
    )
    parser.add_argument("--cross-compile",
      action="store_true", dest='cross_compile', default=False,
      help="Generate bindings and target Windows cross platform build, this will result in different Size/SSize types. This also implies skipping the binding compilation.",
    )

    parser.add_argument("-j", "--jobs",
      default=1,
      type=int,
      help="Number of processors to use on when building. (default: %(default)s)",
    )

    parser.add_argument('-v', "--verbose",
      action="store_true", default=False,
      help="Generate verbose output.",
    )
    
    options = parser.parse_args(args=args[1:])
    logging.info("Options:\n%r", options)

    if options.verbose:
        logging.root.setLevel(logging.DEBUG)

    if options.target and options.one:
        raise ValueError("Can not specify both --one and --target.")

    #Expand relative to absolute bindings target path
    package_path = os.path.abspath(options.package_path)
    bindings_path = path.join(package_path, "rosetta")
    bindings_header_path = path.join(package_path, "rosetta/include")

    # Resolve working directories from repository base directory
    repository_base_dir = subprocess.check_output('git rev-parse --show-toplevel', shell=True).strip()
    logging.info("Resolved base dir: %s", repository_base_dir)

    rosetta_source_path = path.join(repository_base_dir, "source")
    logging.info("Resolved source dir: %s", rosetta_source_path)

    #Switch to python build directory
    os.chdir(script_root_dir)

    # Resolve absolute path of the bindings 'src' directory
    bindings_src_path = os.path.abspath(path.join(script_root_dir, "src"))


    if not os.path.isdir(bindings_path):
        os.makedirs(bindings_path)

    if not os.path.isdir(bindings_header_path):
        os.makedirs(bindings_header_path)

    source_file_patterns = ["*.py"]
    execute(
        'Copy init script and python files...',
        'cp %s %s/' %
            (" ".join(os.path.join("src", p) for p in source_file_patterns), bindings_path),
    )

    prepareBoostLibs(options, bindings_path)
    preparePythonLibs(options, bindings_path)

    rosetta_libs = prepareRosettaLibs(options, rosetta_source_path, bindings_path)

    os.chdir(path.join(rosetta_source_path, "src"))

    if options.one:
        build_targets = options.one
    elif options.target:
        build_targets = candidate_target_modules[:candidate_target_modules.index(options.target)+1]
    else:
        build_targets = candidate_target_modules

    logging.info('Building namespaces: %s', build_targets)
    stageHeaders( options, build_targets, bindings_header_path )
    buildModules(options, build_targets, bindings_src_path, bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, rosetta_libs=rosetta_libs)
    stageStaticFiles(repository_base_dir, script_root_dir, package_path)

def execute(message, command_line, return_=False):
    logging.info("%s:\n%s", message, command_line)
    po = subprocess.Popen(command_line, bufsize=0, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

    output = po.stdout.read()

    logging.debug("Output:\n%s\n", command_line, output)

    if po.returncode is None:
        po.wait()
    res = po.returncode

    if res:
        logging.error("Error: " + command_line)
        logging.error("Error: " + command_line + "\n" + output.decode("ascii", "ignore"))
        if not return_:
            raise ValueError(command_line)

    return res

def getCompilerOptions(options):
    if current_platform != 'macos':
        add_option = '-ffloat-store -ffor-scope'
        if  current_platform_bits == '32':
            add_option += ' -malign-double'
        else:
            add_option += ' -fPIC'

    elif options.compiler == 'clang':
        add_option = '-pipe -ffast-math -funroll-loops -finline-functions -fPIC'
    else:
        add_option = '-pipe -ffor-scope -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s -fPIC'

    #if current_platform == 'cygwin' : add_option =''
    add_option += ' -DBOOST_PYTHON_MAX_ARITY=25 -DPYROSETTA'

    if options.debug_bindings:
        add_option += ' -g -DDEBUG -O0'
    else:
        add_option += ' -DNDEBUG -O3'

    if not options.numpy_support:
        add_option += ' -DPYROSETTA_NO_NUMPY'

    if options.compiler == 'clang':
        add_option += ' -w'

    return add_option

def getLinkerOptions(options):
    ''' Return appropriate linking options based on platform info
    '''
    add_loption = ''
    #if current_platform == 'linux':
    if current_platform != 'macos':  # Linux and cygwin...
        add_loption += '-shared'
        #if current_platform_bits == '32' and current_platform != 'cygwin': add_loption += ' -malign-double'
        if current_platform_bits == '32' : add_loption += ' -malign-double'
    else: add_loption = '-dynamiclib -Xlinker -headerpad_max_install_names'

    return add_loption

def getPlatformIncludePath(options):
    if current_platform == "macos":
        return '../src/platform/macos'
    else:
        return '../src/platform/linux'

def stageHeaders(options, module_paths, dest_inc_path):
    """Recursively stage all headers for given modules."""

    dir_list = []
    for module_path in module_paths:
        for dir_name, _, _ in os.walk(module_path):
            if exclude.isBanned(dir_name):
                logging.info('Skipping banned directory %s.', dir_name)
                continue

            dir_list.append((dir_name, path.join(dest_inc_path, dir_name) ))

    if not options.jobs or options.jobs <= 1:
        map( copy_header_directory, dir_list )
        copy_header_directory((path.join(getPlatformIncludePath(options), "platform"), path.join(dest_inc_path, "platform")))
    else:
        work_pool = multiprocessing.Pool( processes = options.jobs, initializer=init_worker )
        try:
            work_pool.map( copy_header_directory, dir_list )
            work_pool.apply(copy_header_directory, ((path.join(getPlatformIncludePath(options), "platform"), path.join(dest_inc_path, "platform")),))
        except KeyboardInterrupt:
            work_pool.terminate()
            work_pool.join()
            raise
        finally:
            work_pool.close()
            work_pool.join()

def copy_header_directory((source_header_dir, target_header_dir)):
    """Copy all namespace headers into output include directory."""

    header_glob = path.join(source_header_dir, "*.hh")
    all_headers = [p for p in glob.glob(header_glob) if path.isfile(p)]
    logging.debug("Copying headers: %s\n%s", header_glob, pprint.pformat(all_headers) )

    if all_headers and not path.exists(target_header_dir):
        os.makedirs(target_header_dir)

    for header_file in all_headers:
        target_header_file = path.join(target_header_dir, path.basename(header_file))
        if not os.path.isfile(target_header_file) or os.path.getmtime(target_header_file) < os.path.getmtime(header_file):
            logging.debug("Copying header: %s", header_file, )
            shutil.copyfile( header_file, target_header_file )

def buildModules(options, module_paths, py_src_path, dest, include_paths, libpaths, runtime_libpaths, rosetta_libs):
    ''' recursive build buinding for given dir name, and store them in dest.'''

    dir_list = []
    for module_path in module_paths:
        for dir_name, _, files in os.walk(module_path):
            if exclude.isBanned(dir_name):
                logging.info('Skipping banned directory %s.', dir_name)
                continue

            dir_list.append( (dir_name, files) )

    # sort dirs by number of files, most populated first. This should improve speed of multi-thread builds
    dir_list.sort(key=lambda x: -len(x[1]))

    logging.info("Building directories:\n%s", pprint.pformat(dir_list, width=80))

    module_builders = []
    for dir_name, _ in dir_list:
        module_builders.append(ModuleBuilder(options, dir_name, py_src_path, dest, include_paths, libpaths, runtime_libpaths, rosetta_libs))

    if not options.jobs or options.jobs <= 1:
      generate_results = map( perform_build, module_builders)
    else:
        work_pool = multiprocessing.Pool( processes = options.jobs, initializer=init_worker )
        try:
            for res in work_pool.imap_unordered(perform_build, module_builders):
                pass
        except KeyboardInterrupt:
            work_pool.terminate()
            work_pool.join()
            raise
        except Exception:
            work_pool.terminate()
            work_pool.join()
            raise
        finally:
            work_pool.close()
            work_pool.join()

# Global utility functions for use with multiprocessing Pool
def perform_build(module_builder):
    module_builder.generateBindings()
    module_builder.compileBindings()
    module_builder.linkBindings()

def prepareRosettaLibs(options, rosetta_source_path, bindings_path):
    #mode = 'pyrosetta_debug' if options.debug else 'pyrosetta'

    cmake_path = os.path.join(rosetta_source_path, "cmake")
    build_path = os.path.join(cmake_path, "build_pyrosetta")

    subprocess.check_call(["./make_project.py", "all"], cwd=cmake_path )
    subprocess.check_call(["cmake", "-G", "Ninja"], cwd=build_path)

    # Cleanup linked libraries.
    map( os.remove, glob.glob( path.join(build_path, "*." + current_lib_suffix)))

    # Build libraries
    subprocess.check_call(["ninja", candidate_target_library[options.target if options.target else "protocols"]
], cwd=build_path)

    source_libs = glob.glob( path.join(build_path, "*." + current_lib_suffix))
    for s in source_libs:
        shutil.copy( s, bindings_path )

    # Trim off extension and lib prefix
    return [path.splitext(path.basename(s))[0][3:] for s in source_libs]

def prepareBoostLibs(options, bindings_path):
    """Identify and copy boost library into the bindings path."""

    if options.package_boost:
        if not options.boost_path:
            raise ValueError("Must specify options.boost_path if options.package_boost is set.")

        # Identify and copy boost library into target directory
        search_glob = "%s/lib/lib%s.%s*" % (options.boost_path, options.boost_lib, current_lib_suffix)
        logging.debug("prepareBoostLibs searching: %s", search_glob)
        boost_libs = glob.glob(search_glob)
        logging.info("prepareBoostLibs found boost lib: %s", boost_libs)

        if len(boost_libs) == 0:
            raise ValueError("No valid boost library found: %s" % search_glob)

        for boost_lib in boost_libs:
            execute("Copying boost lib: %s" % boost_lib, "cp -P %s %s" % (boost_lib, bindings_path) )
    else:
        if options.boost_path:
            # Append boost include path to options.L
            boost_lib_path = os.path.join(options.boost_path, "lib")
            logging.info("prepareBoostLibs using boost lib path: %s", boost_lib_path)
            options.L.append(boost_lib_path)

    if options.boost_path:
        # Append boost include path to options.I
        boost_include_path = os.path.join(options.boost_path, "include")
        logging.info("prepareBoostLibs using boost include path: %s", boost_include_path)
        options.I.append(boost_include_path)

def preparePythonLibs(options, bindings_path):
    """Identify python include paths.

    Must specify lib and include paths for target python interpreter, as well as numpy include path."""

    python_include_path = None
    python_lib_path = None

    if options.python_path:
        python_include_path = os.path.join(options.python_path, "include", options.python_lib)
        if not os.path.exists(python_include_path):
            logging.warning("Invalid python_path & python_version, include path does not exist: %s", python_include_path)
            logging.warning("Falling back to distutils.sysconfig.get_python_inc.")
            python_include_path = None

        python_lib_path = os.path.join(options.python_path, "lib", options.python_lib)
        if not os.path.exists(python_include_path):
            logging.warning("Invalid python_path & python_version, lib path does not exist: %s", python_lib_path)
            logging.warning("Falling back to distutils.sysconfig.get_python_inc.")
            python_lib_path = None

    import distutils.sysconfig
    if not python_include_path:
        python_include_path = distutils.sysconfig.get_config_var("INCLUDEPY")
    if not python_lib_path:
        python_lib_path = distutils.sysconfig.get_config_var("LIBDIR")

    python_numpy_include_path = None

    if options.python_path:
        python_numpy_include_path = os.path.join(
            options.python_path, "lib", options.python_lib,
            "site-packages", "numpy", "core", "include")

        if not os.path.exists(python_numpy_include_path):
            logging.warning("Numpy not installed under target python_path: %s", python_numpy_include_path)
            logging.warning("Falling back to numpy.get_include.")
            python_numpy_include_path = None

    if not python_numpy_include_path:
        import numpy
        python_numpy_include_path = numpy.get_include()

    logging.info("Using python include path: %s", python_include_path)
    options.I.append(python_include_path)
    logging.info("Using python lib path: %s", python_lib_path)
    options.L.append(python_lib_path)

    logging.info("Using python numpy include path: %s", python_numpy_include_path)
    options.I.append(python_numpy_include_path)

def copy_tree_contents(src, dst, symlinks=False, ignore=None):
    """Recursively copy contents of src into dst."""
    from shutil import Error, WindowsError, copy2, copystat

    names = os.listdir(src)
    if ignore is not None:
        ignored_names = ignore(src, names)
    else:
        ignored_names = set()

    if not path.exists(dst):
        os.makedirs(dst)

    errors = []
    for name in names:
        if name in ignored_names:
            continue
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copy_tree_contents(srcname, dstname, symlinks, ignore)
            else:
                copy2(srcname, dstname)
        except (IOError, os.error) as why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except Error as err:
            errors.extend(err.args[0])
    try:
        copystat(src, dst)
    except WindowsError:
        # can't copy file access times on Windows
        pass
    except OSError as why:
        errors.extend((src, dst, str(why)))
    if errors:
        raise Error(errors)

def stageStaticFiles(repository_root_dir, script_root_dir, target_dir):
    logging.info("stageStaticFiles: %s", locals())

    copy_tree_contents( path.join(script_root_dir, "static"), target_dir)

    if not path.lexists( path.join(target_dir, "rosetta/database") ):
        os.symlink(
            path.relpath( path.join(repository_root_dir, "database"), path.join(target_dir, "rosetta") ),
            path.join(target_dir, "rosetta/database"))

class ModuleBuilder:
    def __init__(self, options, namespace_path, py_src_path, dest, include_paths, libpaths, runtime_libpaths, rosetta_libs):
        ''' Non recursive build buinding for given dir name, and store them in dest.
            options - ModuleBuilder global options object
            namespace_path - relative path to namespace
            py_src_path - path to bindings src path, containing byhand and .py source files
            dest - path to root file destination, actual dest will be dest + path

            return dict of a newly found headers.

            This is a class because we want to generate path/name only once etc.
        '''
        logging.debug( "ModuleBuilder.init(%s)", locals())

        self.options = options
        self.namespace_path = namespace_path
        self.py_src_path  = py_src_path
        self.dest = dest

        # Creating list of headers
        self.headers = [p for p in glob.glob(path.join(self.namespace_path, "*.hh")) if path.isfile(p) and not p.endswith("fwd.hh")]
        self.headers.sort()
        for h in self.headers[:]:
            if exclude.isBanned(h):
                logging.info("Skipping banned header: %s", h)
                self.headers.remove(h)

        logging.debug("Building headers:\n%s", pprint.pformat(self.headers))

        self.fname_base = path.join(self.dest, self.namespace_path)
        if not os.path.isdir(self.fname_base):
            os.makedirs(self.fname_base)

        #Create blank __init__.py for module
        with open( path.join( self.fname_base, "__init__.py"), "w") as f:
            pass

        if not self.headers:
            return

        # Resolve platform specific include path, used in compile & gccxml passes
        self.platform_include_path = getPlatformIncludePath(self.options)

        # Resolve module include paths
        self.include_paths = list(include_paths)
        self.include_paths.append( self.py_src_path )
        #Building from source/src
        self.include_paths.append( "../src" )
        self.include_paths.append(self.platform_include_path)

        self.libpaths = [dest] + libpaths

        # Resolve relative path from namespace outpt dir to binding root directory.
        self.dest_origin_rpath = path.relpath(".",self.namespace_path)
        self.runtime_libpaths = runtime_libpaths +  ["'%s'" % path.join("$ORIGIN", self.dest_origin_rpath) ]

        self.rosetta_libs = rosetta_libs

        self.cpp_defines = '-DPYROSETTA -DBOOST_SYTEM_ -DBOOST_NO_MT -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DPYROSETTA_DISABLE_LCAST_COMPILE_TIME_CHECK -DPTR_BOOST -DPTR_MODERN -DUNUSUAL_ALLOCATOR_DECLARATION'

        # See http://lists.mech.kuleuven.be/pipermail/orocos-dev/2014-April/012832.html
        self.gccxml_options = '-DBOOST_THREAD_DONT_USE_CHRONO -DEIGEN_DONT_VECTORIZE'
        if current_platform == 'macos':
            self.gccxml_options += ' --gccxml-compiler g++-4.2'
            self.gccxml_options += ' -march=nocona'
            self.gccxml_options += ' -fpermissive'

        self.cc_files = []
        self.add_option  = getCompilerOptions(options)
        self.add_loption = getLinkerOptions(options)

        self.by_hand_beginning_file = path.join(self.py_src_path, self.namespace_path,
                                                    '_%s__by_hand_beginning.cc' % self.namespace_path.split('/')[-1])
        if os.path.isfile(self.by_hand_beginning_file):
            logging.info("Loading by_hand beginning file: %s", self.by_hand_beginning_file)
            self.by_hand_beginning = open(self.by_hand_beginning_file).read()
        else:
            self.by_hand_beginning = ""

        self.by_hand_ending_file = path.join(self.py_src_path, self.namespace_path,
                                                    '_%s__by_hand_ending.cc' % self.namespace_path.split('/')[-1])

        if os.path.isfile(self.by_hand_ending_file):
            logging.info("Loading by_hand ending file: %s", self.by_hand_ending_file)
            self.by_hand_ending = open(self.by_hand_ending_file).read()
        else:
            self.by_hand_ending = ""

        self.all_at_once_base = '__' + self.namespace_path.split('/')[-1] + '_all_at_once_'
        self.all_at_once_source_cpp = self.fname_base + '/' + self.all_at_once_base + '.source.cc'
        self.all_at_once_cpp = self.fname_base + '/' + self.all_at_once_base + '.'
        self.all_at_once_obj = self.fname_base + '/' + self.all_at_once_base + '.'
        self.all_at_once_xml = self.fname_base + '/' + self.all_at_once_base + '.xml'
        self.all_at_once_lib = self.fname_base + '/' + self.all_at_once_base + '.so'
        self.all_at_once_json = self.fname_base + '/' + self.all_at_once_base + '.json'
        self.all_at_once_relative_files = []

    def generateBindings(self):
        ''' This function only generate XML file, parse it and generate list of sources that saved in sources.json. '''

        # Setup __init__ for module
        src_init_file = path.join( self.py_src_path, self.namespace_path , '__init__.py')
        dest_init_file = path.join( self.dest, self.namespace_path , '__init__.py')

        if os.path.isfile(src_init_file):
            with open(src_init_file) as i, open(dest_init_file, "a") as o:
                o.write(i.read())
        elif getattr( self, "all_at_once_base", None ):
            with open( dest_init_file, "a") as f:
                f.write("from %s import *\n" % self.all_at_once_base)

        if not self.headers:
            return

        xml_recompile = False

        if not self.options.update:
            xml_recompile = True

        if not xml_recompile:
            try:
                if os.path.isfile(self.by_hand_beginning_file) and os.path.getmtime(self.by_hand_beginning_file) > os.path.getmtime(self.all_at_once_json):
                    xml_recompile = True
                elif os.path.isfile(self.by_hand_ending_file) and os.path.getmtime(self.by_hand_ending_file) > os.path.getmtime(self.all_at_once_json):
                    xml_recompile = True

                for header_file in self.headers:
                    if os.path.getmtime(header_file) > os.path.getmtime(self.all_at_once_json):
                        xml_recompile = True
                    else:
                        pass

                if not os.path.exists(self.all_at_once_cpp+'0.cpp'):
                    xml_recompile = True

            except os.error:
                xml_recompile = True

        for header_file in self.headers:
            hbase = header_file.split('/')[-1][:-3]
            hbase = hbase.replace('.', '_')

            fname = path.join(self.fname_base, '_' + hbase + '.cc')
            self.cc_files.append(fname)

            source_fwd_hh = header_file.replace('.hh', '.fwd.hh')
            source_hh = header_file
            source_cc = header_file.replace('.hh', '.cc')

            self.all_at_once_relative_files.extend( [source_fwd_hh, source_hh, source_cc] )  # just collecting file names...

        with open(self.all_at_once_source_cpp, 'w') as f:
            for header_file in self.headers:
                f.write('#include <%s>\n' % header_file)

        if not xml_recompile:
            logging.debug("Skipping generate pass: %s", self.namespace_path)
            return

        if os.path.isfile(self.all_at_once_lib):
            os.remove(self.all_at_once_lib)

        gen_xml_result = execute('Generating XML representation...',
            "gccxml "
            "-fxml=%(out_xml)s "
            "%(src_cpp)s "
            "%(cpp_defines)s "
            "-I. -I../external/include -I../external/boost_1_55_0  -I../external/dbio "
            "-I%(platform_include_path)s "
            "%(extra_include_paths)s "
            "%(gccxml_options)s "
            "-DBOOST_NO_INITIALIZER_LISTS" % dict(
                gccxml_options =self.gccxml_options,
                src_cpp = self.all_at_once_source_cpp,
                out_xml = self.all_at_once_xml,
                cpp_defines = self.cpp_defines,
                platform_include_path = self.platform_include_path,
                extra_include_paths = " ".join("-I%s" % p for p in self.include_paths)),
            self.options.continue_on_error)

        if gen_xml_result: 
            return

        namespaces_to_wrap = ['::'+self.namespace_path.replace('/', '::')+'::']

        code = tools.CppParser.parseAndWrapModule(
                                self.all_at_once_base,
                                namespaces_to_wrap,
                                self.all_at_once_xml,
                                self.all_at_once_relative_files,
                                max_funcion_size=self.options.max_function_size,
                                by_hand_beginning=self.by_hand_beginning,
                                by_hand_ending=self.by_hand_ending)

        logging.info('Getting include list...')
        includes = exclude.getIncludes(self.headers)

        logging.info('Finalizing[%s]', len(code))
        source_list = []
        for i in range( len(code) ):
            all_at_once_N_cpp = self.all_at_once_cpp+'%s.cpp' % i
            all_at_once_N_obj = self.all_at_once_obj+'%s.o' % i
            source_list.append((all_at_once_N_cpp, all_at_once_N_obj))

            if os.path.isfile(all_at_once_N_obj):
                os.remove(all_at_once_N_obj)

            for fl in self.headers:
                code[i] = '#include <%s>\n' % fl + code[i]

            with open(all_at_once_N_cpp, 'w') as f:
                f.write(code[i])

            exclude.finalize2(all_at_once_N_cpp, self.dest, self.namespace_path, module_name=self.all_at_once_base, includes=includes)
        logging.info('Done!')

        json.dump(source_list, file(self.all_at_once_json, 'w') )

    def compileBindings(self):
        ''' Build early generated bindings.
        '''
        if not self.headers:
            return
        source_list = json.load( file(self.all_at_once_json) )

        recompile = False

        if not self.options.update:
            recompile = True

        if not recompile:
            for (all_at_once_N_cpp, all_at_once_N_obj) in source_list:
                if not os.path.isfile(all_at_once_N_obj)  or  os.path.getmtime(all_at_once_N_cpp) > os.path.getmtime(all_at_once_N_obj):
                    recompile = True
                    break

        if not recompile:
            logging.debug("Skipping compile pass: %s", self.namespace_path)
            return

        for (all_at_once_N_cpp, all_at_once_N_obj) in source_list:
            compiler_cmd = "%(compiler)s %(fname)s -o %(obj_name)s -c %(add_option)s %(cpp_defines)s -I../external/include -I../external/boost_1_55_0 -I../external/dbio %(include_paths)s "
            compiler_dict = dict(
                    add_option=self.add_option,
                    fname=all_at_once_N_cpp,
                    obj_name=all_at_once_N_obj,
                    include_paths=" ".join(["-I%s" % p for p in self.include_paths]),
                    compiler=self.options.compiler, cpp_defines=self.cpp_defines)

            execute("Compiling...", compiler_cmd % compiler_dict, self.options.continue_on_error)

    def linkBindings(self):
        ''' Build early generated bindings.
        '''
        if not self.headers:
            return

        source_list = json.load( file(self.all_at_once_json) )

        relink = False
        if not self.options.update:
            relink = True

        if not relink:
            if not os.path.isfile( self.all_at_once_lib ):
                relink = True
            else:
                for (all_at_once_N_cpp, all_at_once_N_obj) in source_list:
                    if os.path.getmtime( all_at_once_N_obj) > os.path.getmtime(self.all_at_once_lib):
                        relink = True
                        break

        if not relink:
            logging.debug("Skipping link pass: %s", self.namespace_path)
            return

        objs_list = map(lambda x:x[1], source_list)
        linker_cmd = "cd %(dest)s/../ && %(compiler)s %(obj)s %(add_option)s %(rosetta_libs)s -lstdc++ -lz -l%(python_lib)s \
                        -l%(boost_lib)s %(libpaths)s -Wl,%(runtime_libpaths)s -o %(dst)s"
        linker_dict = dict(
                add_option=self.add_loption,
                obj=' '.join(objs_list),
                dst=self.all_at_once_lib,
                libpaths=' '.join(["-L%s" % p for p in self.libpaths]),
                runtime_libpaths=','.join(['-rpath,%s' % p for p in self.runtime_libpaths]),
                rosetta_libs = " ".join("-l%s" % l for l in self.rosetta_libs),
                dest=self.dest,
                boost_lib=self.options.boost_lib,
                python_lib=self.options.python_lib,
                compiler=self.options.compiler)

        execute("Linking...", linker_cmd % linker_dict, self.options.continue_on_error)

if __name__ == "__main__":
    main(sys.argv)
