#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true: :collapseFolds=10:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   BuildBuindings.py
## @brief  Build Python buidings for mini
## @author Sergey Lyskov

import logging
logging.basicConfig(level=logging.DEBUG, format="%(levelname)s %(message)s")
logger = logging.getLogger("BuildPackagedBindings")

import os, re, sys, time, commands, shutil, platform, os.path, gc, json, glob
from os import path
import subprocess #, errno

# Expand path to include tools subdir
script_root_dir = path.dirname(path.realpath(__file__))
logging.info("Script root directory: %s", script_root_dir)
sys.path.append( path.dirname( script_root_dir  ))

# Create global 'Platform' that will hold info of current system
if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
elif sys.platform == "win32" : Platform = "windows"
else: Platform = "_unknown_"
PlatformBits = platform.architecture()[0][:2]

import exclude

import tools.DoxygenExtractorPyPP
import tools.CppParser

from argparse import ArgumentParser

from collections import namedtuple
JobTuple = namedtuple("JobTuple",["pid", "tag"])

Jobs = []  # Global list of NameTuples  (pid, tag)

candidate_target_modules = ['utility', 'numeric', 'basic', 'core', 'protocols']

def main(args):
    ''' Script to build mini Python buidings.
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
        help="Target build namespace, dependent namespaces built.")

    parser.add_argument("--BuildMiniLibs",
      default=True,
      action="store_true",
      help="Build rosetta libraries before generating bindings."
      )

    parser.add_argument("-d",
      action="store_false", dest="BuildMiniLibs",
      help="Disable building of mini libs.",
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
      help="Build bindings with -DDEBUG.",
      )

    parser.add_argument("--numpy_support",
      action="store_true", default=True,
      help="Enable numpy type conversion support.",
      )

    parser.add_argument("--no-numpy_support",
      action="store_true", default=False, dest="numpy_support",
      help="Disable numpy type conversion support.",
      )

    parser.add_argument("--package_path",
        action="store", default="pyrosetta",
        help="Target package directory for bindings. (default: %(default)s)",)

    parser.add_argument("--bindings_path",
      action="store", default=None,
      help="Output directory for bindings. (default: %(default)s)",
      )

    parser.add_argument("--continue",
      default=False,
      action="store_true", dest="continue_",
      help="Debug only. Continue building after encounter the error.",
      )

    parser.add_argument("--all", default=True,
      action="store_true", dest="build_all",
      help="Experimental. Build bindings for all source avaliable.",
      )

    parser.add_argument("--gccxml",
      default='gccxml',
      action="store",
      help="Path to gccxml executable. (default: %(default)s)",
      )

    parser.add_argument("--compiler",
      default='gcc',
      action="store",
      help="Default compiler that will be used to build PyRosetta. (default: %(default)s)",
      )

    parser.add_argument("--gccxml-compiler",
      default='',
      action="store",
      help="Default compiler that will be used in GCCXML. Default is empty string which usually imply 'gcc'.",
      )

    parser.add_argument(
            "--boost_path",
            action="store",
            required=True,
            help="Path to boost install prefix.")

    parser.add_argument("--boost_lib",
      default='boost_python',
      action="store",
      help="Name of boost dynamic library. (default: %(default)s)",
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
        required=True,
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

    parser.add_argument("-p", "--parsing-jobs",
      default=1,
      type=int,
      help="Number of processors to use for parsing when building. WARNING: Some namespace will consume huge amount of memory when parsing (up to ~4Gb), use this option with caution! (default: %(default)s)",
    )

    parser.add_argument('--color',
      action="store_true", default=False,
      help="Color output.",
    )

    parser.add_argument('--no-color',
      action="store_false", dest="color", default=False,
      help="Disable color output.",
    )

    parser.add_argument('-v', "--verbose",
      action="store_true", default=False,
      help="Generate verbose output.",
    )

    options = parser.parse_args(args=args[1:])
    logger.info("Options:\n%r", options)

    global Options
    Options = options

    if Options.parsing_jobs > Options.jobs:
        Options.parsing_jobs = Options.jobs  # Seriously now...

    if Options.target and Options.one:
        raise ValueError("Can not specify both --one and --target.")

    #Expand relative to absolute bindings target path
    package_path = os.path.abspath(options.package_path)
    bindings_path = path.join(package_path, "rosetta")

    # Resolve working directories from repository base directory
    repository_base_dir = subprocess.check_output('git rev-parse --show-toplevel', shell=True).strip()
    logger.info("Resolved base dir: %s", repository_base_dir)

    mini_path = path.join(repository_base_dir, "source")
    logging.info("Resolved source dir: %s", mini_path)

    #Switch to python build directory
    os.chdir(script_root_dir)

    # Resolve absolute path of the bindings 'src' directory
    bindings_src_path = os.path.abspath(path.join(script_root_dir, "src"))


    if not os.path.isdir(bindings_path):
        os.makedirs(bindings_path)

    source_file_patterns = ["*.py"]
    execute(
            'Copy init script and python files...', 'cp %s %s/' %
            (" ".join(os.path.join("src", p) for p in source_file_patterns), bindings_path),
            verbose=Options.verbose)

    prepareBoostLibs(bindings_path)
    preparePythonLibs(bindings_path)

    if options.BuildMiniLibs:
        prepareMiniLibs(mini_path, bindings_path)

    os.chdir(path.join(mini_path, "src"))

    if options.one:
        build_targets = options.one
    elif options.target:
        build_targets = candidate_target_modules[:candidate_target_modules.index(options.target)+1]
    else:
        build_targets = candidate_target_modules

    logger.info('Building namespaces: %s', build_targets)
    buildModules(build_targets, bindings_src_path, bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, gccxml_path=options.gccxml)

    error = False
    for j in Jobs:
        try:
            r = os.waitpid(j.pid, 0)  # waiting for all child process to termintate...
            if r[1] :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                error = True

        except OSError:
            error = True

    if error:
        logging.error('Some of the build scripts return an error, PyRosetta build failed!')
        sys.exit(1)

    stageStaticFiles(repository_base_dir, script_root_dir, package_path)

    logging.info("Done!")

def execute(message, command_line, return_=False, untilSuccesses=False, print_output=True, verbose=True):
    if verbose:
        print message
        print command_line

    while True:
        po = subprocess.Popen(command_line, bufsize=0, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

        f = po.stdout

        output = ''
        for line in f:
            if print_output:
                print line,
            output += line
            sys.stdout.flush()
        f.close()

        if po.returncode is None:
            po.wait()
        res = po.returncode

        if res and untilSuccesses:
            print "Error while executing %s: %s\n" % (message, output)
            print "Sleeping 60s... then I will retry..."
            time.sleep(60)
        else:
            break

    if res:
        if print_output:
            print "\nEncounter error while executing: " + command_line
        if not return_:
            sys.exit(1)

    if return_ == 'output':
        return output
    else:
        return res

def print_(msg, color=None, background=None, bright=False, blink=False, action='print', endline=True):
    ''' print string with color and background. Avoid printing and return results str instead if action is 'return'. Also check for 'Options.no_color'
    '''
    colors = dict(black=0, red=1, green=2, yellow=3, blue=4, magenta=5, cyan=6, white=7)  # standard ASCII colors

    if 'Options' in globals()  and  hasattr(Options, 'color')  and  not Options.color: s = str(msg)
    else:
        s  = ['3%s' % colors[color] ] if color else []
        s += ['4%s' % colors[background] ] if background else []
        s += ['1'] if bright else []
        s += ['5'] if blink else []

        s = '\033[%sm%s\033[0m%s' % (';'.join(s), msg, '\n' if endline else '')

    if action == 'print': sys.stdout.write(s)
    else: return s

def mFork(tag=None, overhead=0):
    ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
    '''
    #print_('Groups:%s' % os.getgroups(), color='cyan')
    while len(Jobs) >= Options.jobs + overhead:
        for j in Jobs[:] :
            r = os.waitpid(j.pid, os.WNOHANG)
            if r == (j.pid, 0):  # process have ended without error
                Jobs.remove(j)
            elif r[0] == j.pid :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                for j in Jobs:
                    try:
                        os.waitpid(j.pid, 0)
                    except OSError: pass

                print_('Some of the build scripts return an error, PyRosetta build failed!', color='red', bright=True)
                sys.exit(1)

        if len(Jobs) >= Options.jobs + overhead: time.sleep(.2)

    sys.stdout.flush()
    pid = os.fork()
    if pid: Jobs.append( JobTuple(pid, tag) )  # We are parent!
    return pid

def mWait(tag=None, all_=False):
    ''' Wait for process tagged with 'tag' for completion
    '''
    while True :
        #print 'Waiting for %s: ' % tag, Jobs
        for j in [ x for x in Jobs if x.tag==tag or all_==True]:
            #print 'Waiting2: ', Jobs
            #try:
            r = os.waitpid(j.pid, os.WNOHANG)
            if r == (j.pid, 0):  # process have ended without error
                Jobs.remove(j)
            elif r[0] == j.pid :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                for j in Jobs:
                    try:
                        os.waitpid(j.pid, 0)
                    except OSError: pass

                print_('Some of the build scripts return an error, PyRosetta build failed!', color='red', bright=True)
                sys.exit(1)
            else: time.sleep(.2);  break
            '''
            except OSError, e:
                if e.errno == errno.ESRCH:  # process already got closed, we assume that this is done by child process and any error will be reported by proc. that closed it
                    Jobs.remove(j) '''

        else: return

def getCompilerOptions():
    #if Platform == 'linux':
    if Platform != 'macos':
        add_option = '-ffloat-store -ffor-scope'
        if  PlatformBits == '32':
            add_option += ' -malign-double'
        else:
            add_option += ' -fPIC'
    elif Options.compiler == 'clang':
        add_option = '-pipe -O3 -ffast-math -funroll-loops -finline-functions -fPIC'
    else:
        add_option = '-pipe -ffor-scope -O3 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s -fPIC'
    #if Platform == 'cygwin' : add_option =''
    add_option += ' -DBOOST_PYTHON_MAX_ARITY=25 -DPYROSETTA'
    add_option += (' -DDEBUG' if Options.debug_bindings else ' -DNDEBUG')

    if not Options.numpy_support:
        add_option += ' -DPYROSETTA_NO_NUMPY'

    if Options.compiler == 'clang':
        add_option += ' -w'

    return add_option

def getLinkerOptions():
    ''' Return appropriate linking options based on platform info
    '''
    add_loption = ''
    #if Platform == 'linux':
    if Platform != 'macos':  # Linux and cygwin...
        add_loption += '-shared'
        #if PlatformBits == '32' and Platform != 'cygwin': add_loption += ' -malign-double'
        if PlatformBits == '32' : add_loption += ' -malign-double'
    else: add_loption = '-dynamiclib -Xlinker -headerpad_max_install_names'

    return add_loption

def buildModules(module_paths, py_src_path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
    ''' recursive build buinding for given dir name, and store them in dest.'''

    dir_list = []
    for module_path in module_paths:
        for dir_name, _, files in os.walk(module_path):
            if exclude.isBanned(dir_name):
                logger.info('Skipping banned directory %s.', dir_name)
                continue

            dir_list.append( (dir_name, files) )

    # sort dirs by number of files, most populated first. This should improve speed of multi-thread builds
    dir_list.sort(key=lambda x: -len(x[1]))

    logging.info("Building directories: %s", dir_list)

    mb = []
    for dir_name, _ in dir_list:
        #print "buildModules(...): '%s', " % dir_name
        dname = dest + '/' + dir_name
        if not os.path.isdir(dname): os.makedirs(dname)

        mb.append( ModuleBuilder(dir_name, py_src_path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path) )
        mb[-1].generateBindings()
        gc.collect()

    mWait(all_=True)  # waiting for all jobs to finish before movinf in to next phase

    for b in mb:
        b.compileBindings()
        gc.collect()

    mWait(all_=True)  # waiting for all jobs to finish before movinf in to next phase

    for b in mb:
        b.linkBindings()
        gc.collect()

def prepareMiniLibs(mini_path, bindings_path):
    mode = 'pyrosetta_debug' if Options.debug else 'pyrosetta'

    if Platform == "macos" and PlatformBits=='32':
        execute("Building mini libraries...", "cd %s && ./scons.py mode=%s arch=x86 arch_size=32 -j%s" % (mini_path, mode, Options.jobs) )
    elif Platform == "macos" and PlatformBits=='64':
        execute("Building mini libraries...", "cd %s && ./scons.py mode=%s -j%s" % (mini_path, mode, Options.jobs) )
    elif Platform == "cygwin":
        execute("Building mini libraries...", "cd %s && ./scons.py mode=%s bin -j%s" % (mini_path, mode, Options.jobs) )
    else:
        execute("Building mini libraries...", "cd %s && ./scons.py mode=%s -j%s" % (mini_path, mode, Options.jobs) )

    # fix this for diferent platform
    if Platform == "linux":
        lib_path = os.path.join('build/src/', mode, 'linux/' , platform.release()[:3], PlatformBits , 'x86/gcc/')
    elif Platform == "cygwin":
        lib_path = os.path.join('build/src/', mode, 'cygwin/1.7/32/x86/gcc/')
    else:
        if Platform == "macos" and PlatformBits=='32':
            lib_path = os.path.join('build/src/', mode, 'macos/10.5/32/x86/gcc/')
        if Platform == "macos" and PlatformBits=='64':
            if platform.release()[:2] == '10':
                lib_path = os.path.join('build/src/', mode,'macos/10.6/64/x86/gcc/')
            elif platform.release()[:2] == '11':
                lib_path = os.path.join('build/src/', mode, 'macos/10.7/64/x86/gcc/')
            else:
                lib_path = os.path.join('build/src/', mode, 'macos/10.8/64/x86/gcc/')

    # now lets add version to lib_path...
    lib_path += execute("Getting GCC version...", 'gcc -dumpversion', return_='output').strip()[0:3] + '/default/'

        #if Platform == "macos" and PlatformBits=='64'  and  platform.release().startswith('11.'): lib_path = 'build/src/pyrosetta/macos/11/64/x86/gcc/'
        #else: lib_path = 'build/src/pyrosetta/macos/10.6/64/x86/gcc/'

    #lib_path += 'static/'
    obj_suffix = '.os'

    # Now the funny part - we link all libs to produce just one lib file...
    all_sources = []
    all_scons_files = [f for f in commands.getoutput('cd ../../ && ls *.src.settings').split() if f not in ['apps.src.settings', 'devel.src.settings', 'pilot_apps.src.settings']]
    for scons_file in all_scons_files:
    #for scons_file in ['ObjexxFCL', 'numeric', 'utility',]:
        f = file('./../../'+scons_file).read();  exec(f)
        for k in sources:
            for f in sources[k]:
                #all_sources.append( scons_file + '/' + k + '/' + f + obj_suffix)
                if not f.endswith('.cu'): all_sources.append( k + '/' + f + obj_suffix)

    #all_sources.remove('protocols/forge/remodel/RemodelDesignMover' + obj_suffix) # <-- I have no idea what gcc does not like this file...

    extra_objs = [  # additioanal source for external libs
        'dbio/cppdb/atomic_counter', "dbio/cppdb/conn_manager", "dbio/cppdb/driver_manager", "dbio/cppdb/frontend",
        "dbio/cppdb/backend", "dbio/cppdb/mutex", "dbio/cppdb/pool", "dbio/cppdb/shared_object", "dbio/cppdb/sqlite3_backend",
        "dbio/cppdb/utils", 'dbio/sqlite3/sqlite3', ]

    all_sources += [ mini_path + '/' + lib_path.replace('/src/', '/external/') + x + obj_suffix for x in extra_objs ]

    objs = ' '.join(all_sources)

    suffix = 'so'
    if Platform == 'cygwin':
        suffix = 'dll'
    if Platform == 'macos':
        suffix = 'dylib'

    mini = os.path.join(bindings_path, 'libmini.' + suffix)

    add_loption = getLinkerOptions()

    execute("Linking mini lib...",
            "cd %(mini_path)s && cd %(lib_path)s && gcc %(add_loption)s \
            %(objs)s -lz -lstdc++ -o %(mini)s" % dict(mini_path=mini_path, lib_path=lib_path, add_loption=add_loption, mini=mini, objs=objs, compiler=Options.compiler)
             )

    if Platform == 'macos':
        #libs = ['libObjexxFCL.dylib', 'libnumeric.dylib', 'libprotocols.dylib', 'libdevel.dylib', 'libutility.dylib', 'libcore.dylib']
        libs = ['libmini.dylib']
        for l in libs:
            execute('Adjustin lib self path in %s' % l, 'install_name_tool -id rosetta/%s %s' % (l, os.path.join(bindings_path, l)) )
            for k in libs:
                execute('Adjustin lib path in %s' % l, 'install_name_tool -change %s rosetta/%s %s' % (os.path.abspath(path.join(mini_path,lib_path,k)), k, os.path.join(bindings_path, l)) )

def prepareBoostLibs(bindings_path):
    """Identify and copy boost library into the bindings path."""

    #TODO alexford factor out suffix gen
    suffix = 'so'
    if Platform == 'cygwin':
        suffix = 'dll'
    if Platform == 'macos':
        suffix = 'dylib'

    # Identify and copy boost library into target directory
    search_glob = "%s/lib/lib%s.%s*" % (Options.boost_path, Options.boost_lib, suffix)
    logging.debug("prepareBoostLibs searching: %s", search_glob)
    boost_libs = glob.glob(search_glob)
    logging.info("prepareBoostLibs found boost lib: %s", boost_libs)

    if len(boost_libs) == 0:
        raise ValueError("No valid boost library found: %s" % search_glob)

    for boost_lib in boost_libs:
        execute("Copying boost lib: %s" % boost_lib, "cp -P %s %s" % (boost_lib, bindings_path) )

    # Append boost include path to Options.I
    boost_include_path = os.path.join(Options.boost_path, "include")
    logging.info("prepareBoostLibs using boost include path: %s", boost_include_path)
    Options.I.append(boost_include_path)

def preparePythonLibs(bindings_path):
    """Identify python include paths."""

    python_include_path = os.path.join(Options.python_path, "include", Options.python_lib)
    if not os.path.exists(python_include_path):
        logger.warning("Invalid python_path & python_version, include path does not exist: %s", python_include_path)
        import distutils
        python_include_path = distutils.sysconfig.get_python_inc()
        logger.warning("Falling back to distutils.sysconfig.get_python_inc: %s", python_include_path)

    python_numpy_include_path = os.path.join(
            Options.python_path, "lib", Options.python_lib,
            "site-packages", "numpy", "core", "include")

    if not os.path.exists(python_numpy_include_path):
        logger.warning("Numpy not installed under target python_path: %s", python_numpy_include_path)
        import numpy
        python_numpy_include_path = numpy.get_include()
        logger.warning("Falling back to numpy.get_include: %s", python_numpy_include_path)

    logger.info("Using python include path: %s", python_include_path)
    Options.I.append(python_include_path)
    logger.info("Using python numpy include path: %s", python_numpy_include_path)
    Options.I.append(python_numpy_include_path)

    python_library_path = os.path.join(Options.python_path, "lib")

    #Just glob on libpythnon<major>.<minor>* to avoid needing to resolve dylib/so/etc...
    shared_lib_glob = os.path.join(python_library_path, "lib" + Options.python_lib + "*")
    if len(glob.glob(shared_lib_glob)) == 0:
        raise ValueError("Unable to resolve python shared library matching include path and version: %s" % shared_lib_glob)

    logger.info("Using python library path: %s", python_library_path)
    Options.L.append(python_library_path)

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
    logger.info("stageStaticFiles: %s", locals())

    copy_tree_contents( path.join(script_root_dir, "static"), target_dir)

    if not path.lexists( path.join(target_dir, "database") ):
        os.symlink(
            path.relpath( path.join(repository_root_dir, "database"), target_dir),
            path.join(target_dir, "database"))

class ModuleBuilder:
    def __init__(self, namespace_path, py_src_path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
        ''' Non recursive build buinding for given dir name, and store them in dest.
            namespace_path - relative path to namespace
            py_src_path - path to bindings src path, containing byhand and .py source files
            dest - path to root file destination, actual dest will be dest + path

            return dict of a newly found headers.

            This is a class because we want to generate path/name only once etc.
        '''
        global Options
        if Options.verbose: print 'CppXML route: ModuleBuilder.init...', path, dest

        self.namespace_path = namespace_path
        self.py_src_path  = py_src_path
        self.dest = dest

        # Creating list of headers
        self.headers = [os.path.join(self.namespace_path, d)
                            for d in os.listdir(self.namespace_path)
                                if os.path.isfile( os.path.join(self.namespace_path, d) )
                                    and d.endswith('.hh')
                                    and not d.endswith('.fwd.hh')
                                    ]

        self.headers.sort()

        for h in self.headers[:]:
            if exclude.isBanned(h):
                if Options.verbose: print "Banning header:", h
                self.headers.remove(h)

        if Options.verbose: print_(self.headers, color='black', bright=True)

        self.fname_base = self.dest + '/' + self.namespace_path

        if not os.path.isdir(self.fname_base): os.makedirs(self.fname_base)

        #print 'Creating __init__.py file...'
        f = file( dest + '/' + self.namespace_path + '/__init__.py', 'w');  f.close()
        if not self.headers:  return   # if source files is empty then __init__.py should be empty too

        def finalize_init(current_fname):
            #print 'Finalizing Creating __init__.py file...'
            f = file( dest + '/' + self.namespace_path + '/__init__.py', 'a');
            f.write('from %s import *\n' % os.path.basename(current_fname)[:-3]);
            f.close()

        self.finalize_init = finalize_init

        # Resolve platform specific include path, used in compile & gccxml passes
        if Platform == "macos":
            self.platform_include_path = '../src/platform/macos'
        else:
            self.platform_include_path = '../src/platform/linux'

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

        self.cpp_defines = '-DPYROSETTA -DBOOST_SYTEM_ -DBOOST_NO_MT -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DPYROSETTA_DISABLE_LCAST_COMPILE_TIME_CHECK'

        self.gccxml_options = ''
        if Options.gccxml_compiler:
            self.gccxml_options += '--gccxml-compiler ' + Options.gccxml_compiler
        elif Platform == 'macos':
            self.gccxml_options += '--gccxml-compiler llvm-g++-4.2'

        if Platform == 'macos':
            self.gccxml_options += ' -march=nocona'

        self.cc_files = []
        self.add_option  = getCompilerOptions()
        self.add_loption = getLinkerOptions()

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
        ''' This function only generate XML file, parse it and generate list of sources that saved in sources.json. We assume that one_lib_file and build_all option is on here.
        '''
        if not self.headers: return

        xml_recompile = False
        for fl in self.headers:
            #print 'Binding:', files
            hbase = fl.split('/')[-1][:-3]
            hbase = hbase.replace('.', '_')
            #print 'hbase = ', hbase
            #if hbase == 'init': hbase='tint'  # for some reason Boost don't like 'init' name ? by hand?

            fname = self.fname_base + '/' + '_' + hbase + '.cc'
            '''
            inc_name =  fname_base + '/' + '_' + hbase + '.hh'
            obj_name =  fname_base + '/' + '_' + hbase + '.o'
            xml_name =  fname_base + '/' + '_' + hbase + '.xml'
            cc_for_xml_name =  fname_base + '/' + '_' + hbase + '.xml.cpp'
            dst_name =  fname_base + '/' + '_' + hbase + '.so'
            if Platform == 'cygwin' : dst_name =  fname_base + '/' + '_' + hbase + '.dll'
            decl_name = fname_base + '/' + '_' + hbase + '.exposed_decl.pypp.txt'
            '''
            self.cc_files.append(fname)

            if Options.update:
                try:
                    if fl == self.headers[0]:  # for first header we additionaly check if 'by_hand' code is up to date
                        if os.path.isfile(self.by_hand_beginning_file) and os.path.getmtime(self.by_hand_beginning_file) > os.path.getmtime(self.all_at_once_json): raise os.error
                        if os.path.isfile(self.by_hand_ending_file) and os.path.getmtime(self.by_hand_ending_file) > os.path.getmtime(self.all_at_once_json): raise os.error

                    if os.path.getmtime(fl) > os.path.getmtime(self.all_at_once_json):
                        xml_recompile = True
                    else:
                        if Options.verbose: print 'File: %s is up to date - skipping' % fl

                except os.error: xml_recompile = True

            if xml_recompile: print_(fl, color='green', bright=True)

            source_fwd_hh = fl.replace('.hh', '.fwd.hh')
            source_hh = fl
            source_cc = fl.replace('.hh', '.cc')

            self.all_at_once_relative_files.extend( [source_fwd_hh, source_hh, source_cc] )  # just collecting file names...


        f = file(self.all_at_once_source_cpp, 'w');
        for fl in self.headers: f.write('#include <%s>\n' % fl);
        f.close()


        #print 'Finalizing Creating __init__.py file...'
        namespace = os.path.basename(self.namespace_path)
        py_init_file = path.join( self.py_src_path, self.namespace_path , '__init__.py')
        if os.path.isfile(py_init_file):
            t = file(py_init_file).read()
        else:
            t = ''

        with open(self.dest + '/' + self.namespace_path + '/__init__.py', 'w') as f:
            f.write(t+'from %s import *\n' % self.all_at_once_base);
            f.close()

        if xml_recompile or (not Options.update):
            start_time = time.time()

            if os.path.isfile(self.all_at_once_lib): os.remove(self.all_at_once_lib)

            def generate():
                if execute('Generating XML representation...', 'gccxml %s %s -fxml=%s %s -I. -I../external/include -I../external/boost_1_55_0  -I../external/dbio -I%s -DBOOST_NO_INITIALIZER_LISTS ' % (self.gccxml_options, self.all_at_once_source_cpp, self.all_at_once_xml, self.cpp_defines, self.platform_include_path), Options.continue_ ): return

                namespaces_to_wrap = ['::'+self.namespace_path.replace('/', '::')+'::']
                # Temporary injecting Mover in to protocols level
                #if path == 'protocols': namespaces_to_wrap.append('::protocols::moves::')

                code = tools.CppParser.parseAndWrapModule(self.all_at_once_base, namespaces_to_wrap, self.all_at_once_xml, self.all_at_once_relative_files, max_funcion_size=Options.max_function_size,
                                                          by_hand_beginning=self.by_hand_beginning, by_hand_ending=self.by_hand_ending)

                print_('Getting include list...', color='black', bright=True)
                includes = exclude.getIncludes(self.headers)

                print_('Finalizing[%s]' % len(code), color='black', bright=True, endline=False);  sys.stdout.flush()
                source_list = []
                for i in range( len(code) ):
                    all_at_once_N_cpp = self.all_at_once_cpp+'%s.cpp' % i
                    all_at_once_N_obj = self.all_at_once_obj+'%s.o' % i
                    source_list.append((all_at_once_N_cpp, all_at_once_N_obj))

                    if os.path.isfile(all_at_once_N_obj): os.remove(all_at_once_N_obj)

                    for fl in self.headers: code[i] = '#include <%s>\n' % fl + code[i]

                    f = file(all_at_once_N_cpp, 'w');  f.write(code[i]);  f.close()

                    exclude.finalize2(all_at_once_N_cpp, self.dest, self.namespace_path, module_name=self.all_at_once_base, add_by_hand = False, includes=includes)
                    print_('.', color='black', bright=True, endline=False); sys.stdout.flush()
                print_(' Done!', color='black', bright=True);

                json.dump(source_list, file(self.all_at_once_json, 'w') )

            if Options.jobs > 1:
                pid = mFork()
                if not pid:  # we are child process
                    generate()
                    sys.exit(0)

            else:
                generate();


    def compileBindings(self):
        ''' Build early generated bindings.
        '''
        if not self.headers: return
        source_list = json.load( file(self.all_at_once_json) )

        recompile = False

        if Options.update:
            for (all_at_once_N_cpp, all_at_once_N_obj) in source_list:
                if not os.path.isfile(all_at_once_N_obj)  or  os.path.getmtime(all_at_once_N_cpp) > os.path.getmtime(all_at_once_N_obj): recompile = True; break

        if recompile or (not Options.update):
            for (all_at_once_N_cpp, all_at_once_N_obj) in source_list: #range( len(code) ):
                start_time = time.time()

                #all_at_once_N_cpp = self.all_at_once_cpp+'%s.cpp' % i
                #all_at_once_N_obj = self.all_at_once_obj+'%s.o' % i

                # -fPIC
                comiler_cmd = "%(compiler)s %(fname)s -o %(obj_name)s -c %(add_option)s %(cpp_defines)s -I../external/include -I../external/boost_1_55_0 -I../external/dbio %(include_paths)s "
                comiler_dict = dict(
                        add_option=self.add_option,
                        fname=all_at_once_N_cpp,
                        obj_name=all_at_once_N_obj,
                        include_paths=" ".join(["-I%s" % p for p in self.include_paths]),
                        compiler=Options.compiler, cpp_defines=self.cpp_defines)

                failed = False

                def compile_():
                    if execute("Compiling...", comiler_cmd % comiler_dict, return_=True):
                        if Options.compiler != 'clang': failed = True
                        elif execute("Compiling...", comiler_cmd % dict(comiler_dict, compiler='gcc'), return_=True): failed = True

                if Options.jobs > 1:
                    pid = mFork(tag=self.namespace_path)
                    if not pid:  # we are child process
                        compile_()
                        sys.exit(0)

                else:
                    compile_()


                if Options.jobs == 1:
                    if Options.continue_ and failed: return new_headers


    def linkBindings(self):
        ''' Build early generated bindings.
        '''
        if not self.headers: return
        source_list = json.load( file(self.all_at_once_json) )

        relink = False

        if Options.update:
            for (all_at_once_N_cpp, all_at_once_N_obj) in source_list:
                if not os.path.isfile(self.all_at_once_lib)  or  os.path.getmtime( all_at_once_N_obj) > os.path.getmtime(self.all_at_once_lib):
                    relink = True
                    break
        else:
            relink = True


        if relink:
            start_time = time.time()

            objs_list = map(lambda x:x[1], source_list)
            linker_cmd = "cd %(dest)s/../ && %(compiler)s %(obj)s %(add_option)s -lmini -lstdc++ -lz -l%(python_lib)s \
                            -l%(boost_lib)s %(libpaths)s -Wl,%(runtime_libpaths)s -o %(dst)s"
            linker_dict = dict(
                    add_option=self.add_loption,
                    obj=' '.join(objs_list),
                    dst=self.all_at_once_lib,
                    libpaths=' '.join(["-L%s" % p for p in self.libpaths]),
                    runtime_libpaths=','.join(['-rpath,%s' % p for p in self.runtime_libpaths]),
                    dest=self.dest,
                    boost_lib=Options.boost_lib,
                    python_lib=Options.python_lib,
                    compiler=Options.compiler)

            def linking():
                if execute("Linking...", linker_cmd % linker_dict, return_= (True if Options.compiler != 'gcc' or Options.continue_ else False) ):
                    if Options.compiler != 'gcc':
                        execute("Linking...", linker_cmd % dict(linker_dict, compiler='gcc'), return_= Options.continue_ )

            if Options.jobs > 1:

                pid = mFork(tag=self.namespace_path+'+linking', overhead=1)  # we most likely can start extra linking process, beceause it depend on compilation to  finish. There is no point of waiting for it...
                if not pid:  # we are child process
                    #mWait(tag=self.path)  # wait for all compilation jobs to finish...
                    linking()
                    sys.exit(0)
            else:
                linking()

if __name__ == "__main__":
    main(sys.argv)
