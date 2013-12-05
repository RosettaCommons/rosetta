#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true: :collapseFolds=10:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   BuildBuindings.py
## @brief  Build Python buidings for Rosetta
## @author Sergey Lyskov

import os, re, sys, time, commands, shutil, platform, os.path, gc, json, types
import subprocess #, errno

# Create global 'Platform' that will hold info of current system
if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
elif sys.platform == "win32" : Platform = "windows"
else: Platform = "_unknown_"
PlatformBits = platform.architecture()[0][:2]


#if Platform != "windows":
#    from pyplusplus import module_builder
#    import pyplusplus

import exclude

import tools.DoxygenExtractorPyPP
import tools.CppParser



from optparse import OptionParser


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            if not i.startswith('__') and not isinstance(getattr(self, i), types.MethodType): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'



# global include dict 'file' --> bool;  True mean include file, False - exclude
# if key is not present that mean that this file/names space is new, it will be excluded and
#    added to IncludeDict.new
IncludeDict = None  # we set it on purpose to None instead of {}, this should cath uninitialize errors

#  Here we will store new files and dirs. This will be saved to 'IncludeDict.new' in the end.
IncludeDictNew = {}

Jobs = []  # Global list of NameTuples  (pid, tag)


def main(args):
    ''' Script to build Rosetta Python buidings.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option('-I',
      default=[],
      action="append",
      help="Additiona path to include (boost, python etc). Specify it only if this libs installed in non standard locations. May specify multiple times.",
    )

    parser.add_option('-L',
      default=[],
      action="append",
      help="Additiona path to libraries (boost, python etc). Specify it only if this libs installed in non standard locations. May specify multiple times.",
    )

    parser.add_option("-1", "--one",
      default=[], action="append",
      help="Build just one namespace instead of whole project, can be specified multiple times.",
    )

    parser.add_option("--BuildMiniLibs",
      default=True,
      )

    parser.add_option("-d",
      action="store_false", dest="BuildMiniLibs",
      help="Disable building of mini libs.",
    )

    parser.add_option("-u", "--update",
      default=False,
      action="store_true",
      help="Debug only. Try to check time stamp of files before building them.",
      )

    parser.add_option("--debug",
      action="store_true", default=False,
      help="Perform a Debug build when possible.",
      )

    parser.add_option("--numpy-support",
      action="store_true", default=False,
      help="Build bindings numpy type conversion support.",
      )

    parser.add_option("--continue",
      default=False,
      action="store_true", dest="continue_",
      help="Debug only. Continue building after encounter the error.",
      )

    parser.add_option("--all", default=True,
      action="store_true", dest="build_all",
      help="Experimental. Build bindings for all source avaliable.",
      )

    parser.add_option("--partial",
      action="store_false", dest="build_all",
      help="Build only limited subset of py-bindings.",
      )

    parser.add_option("--gccxml",
      default='gccxml',
      action="store",
      help="Path to gccxml executable. Default is just 'gccxml'.",
      )

    parser.add_option("--compiler",
      default='gcc',
      action="store",
      help="Default compiler that will be used to build PyRosetta. Default is 'gcc'.",
      )

    parser.add_option("--gccxml-compiler",
      default='',
      action="store",
      help="Default compiler that will be used in GCCXML. Default is empty string which usually imply 'gcc'.",
      )


    parser.add_option("--py-plus-plus", action="store_true", dest="py_plus_plus", default=False,
      help="Use Py++/PyGccXML parser and Python buildes to build PyRosetta [depricated]."
    )

    parser.add_option("--boost-lib",
      default='boost_python-xgcc40-mt-1_38',
      action="store",
      help="Name of boost dynamic library.",
      )

    parser.add_option("--python-lib",
      default='python2.5',
      action="store",
      help="Name of python dynamic library.",
    )

    parser.add_option("--update-IncludeDict",
      action="store_true", dest='sort_IncludeDict', default=False,
      help="Developers only. Save sorter IncludeDict back to file.",
    )

    parser.add_option("--one-lib-file", action="store_true", dest="one_lib_file", default=True,
      help="Generate only one lib file for name spaces [experimental]."
    )

    parser.add_option("--many-lib-files", action="store_false", dest="one_lib_file",
      help="Generate only one lib file for name spaces [experimental]."
    )


    parser.add_option("--one-lib-max-files", default=0, type='int',
      help="Maximum number of files that goes in to one lib file. Default is unlimeted [0]."
    )

    parser.add_option("--max-function-size", default=1024*128, type='int',
      help="Maximum size of function in binding files in bytes."
    )

    parser.add_option("--use-pre-generated-sources",
      default=None,
      action="store",
      help="Mostly for Windows native build: Path to pre-generated PyRosetta C++ source files.",
    )

    parser.add_option("--cross-compile",
      action="store_true", dest='cross_compile', default=False,
      help="Generate bindings and target Windows cross platform build, this will result in different Size/SSize types. This also implies skipping the binding compilation.",
    )

    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="Number of processors to use on when building(default: 1)",
    )

    parser.add_option("-p", "--parsing-jobs",
      default=1,
      type="int",
      help="Number of processors to use for parsing when building. WARNING: Some namespace will consume huge amount of memory when parsing (up to ~4Gb), use this option with caution! (default: 1)",
    )

    parser.add_option("-s", "--sleep",
      default=0,
      type="float",
      help="Sleep number of second between compilation jobs, good for notebook quiet operation! Actual sleep time will be: time_spend_compiling*option_value  (default: 0)",
    )

    parser.add_option('--no-color',
      action="store_false", dest='color', default=True,
      help="Disable color output [Good when piping output to file]",
    )

    parser.add_option('-v', "--verbose",
      action="store_true", default=False,
      help="Generate verbose output.",
    )

    (options, args) = parser.parse_args(args=args[1:])
    global Options;  Options = options

    print "-I", options.I
    print "-L", options.L
    print "-1", options.one
    print "--BuildMiniLibs", options.BuildMiniLibs
    print "--gccxml", options.gccxml
    print "--all", options.build_all
    print "--one-lib-file", options.one_lib_file
    print "--boost-lib", options.boost_lib
    print '--update', options.update
    print '--debug', options.debug
    print '--compiler', options.compiler

    #print '--numpy-support', options.numpy_support

    #print '--cross-compile', options.cross_compile

    if (Options.jobs > 1) and Options.sort_IncludeDict:
        print 'Option --update-IncludeDict could be used only with -j1... exiting...'
        return 1

    if Options.parsing_jobs > Options.jobs: Options.parsing_jobs = Options.jobs  # Seriously now...

    # assuming that we in rosetta/rosetta_source/src/python/bindings directory at that point
    mini_path = os.path.abspath('./../../../')

    bindings_path = os.path.abspath('rosetta')
    #bindings_path += '.debug' if Options.debug else '.release'

    if Options.cross_compile:
        bindings_path = os.path.abspath('rosetta.window')
        if not os.path.isdir(bindings_path): os.makedirs(bindings_path)
        execute('Generating svn_version files...', 'cd ./../../../ && python svn_version.py')  # Now lets generate svn_version.* files and copy it to destination (so windows build could avoid running it).
        shutil.copyfile('./../../core/svn_version.cc', bindings_path + '/svn_version.cc')


    if not os.path.isdir(bindings_path): os.makedirs(bindings_path)
    #shutil.copyfile('src/__init__.py', bindings_path + '/__init__.py')
    execute('Copy init script and additional files...', 'cp src/*.py %s/' % bindings_path, verbose=False)

    if Platform == "windows":  # we dealing with windows native build
        #bindings_path = os.path.abspath('python/bindings/rosetta')
        build_path = os.path.abspath( os.path.join(mini_path, 'build\windows') )
        if Options.debug: build_path    += '_debug'
        BuildRosettaOnWindows(build_path, bindings_path)
        sys.exit(0)

    #if os.path.islink('rosetta') or os.path.isdir('rosetta'): os.remove('rosetta')
    #os.symlink(bindings_path, 'rosetta')

    #bindings_path = os.path.abspath(bindings_path)

    if options.BuildMiniLibs:
        prepareMiniLibs(mini_path, bindings_path)

    os.chdir( './../../' )


    # loadind include file/namespaces dict
    global IncludeDict
    IncludeDict = eval( file('python/bindings/IncludeDict').read() )
    if type( IncludeDict['core'] ) == str:
        for k in IncludeDict: IncludeDict[k] = ({'+':True, '-':False}[IncludeDict[k]], 999, [])

    if options.one:
        print 'Building just following namespaces:', options.one
        for n in options.one:
            #print 'namespace =', n
            buildModule(n, bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, gccxml_path=options.gccxml)

    else:
        buildModules(['utility', 'numeric', 'basic', 'core', 'protocols'],   bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, gccxml_path=options.gccxml)
        '''
        # we want to start with lib that is longest to build - that way we can do multi-core build more efficiently
        buildModules('core',      bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, gccxml_path=options.gccxml)
        buildModules('protocols', bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, gccxml_path=options.gccxml)
        buildModules('utility',   bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, gccxml_path=options.gccxml)
        buildModules('numeric',   bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, gccxml_path=options.gccxml)
        buildModules('basic',     bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L, gccxml_path=options.gccxml)
        '''
    error = False
    for j in Jobs:
        try:
            r = os.waitpid(j.pid, 0)  # waiting for all child process to termintate...
            if r[1] :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                error = True

        except OSError:
            error = True


    if error:
        print_('Some of the build scripts return an error, PyRosetta build failed!', color='red', bright=True)
        sys.exit(1)




    #buildModule(, bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L)
    #buildModule('core/io/pdb', bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L)

    #buildModule('core/kinematics', bindings_path, include_paths=options.I, libpaths=options.L, runtime_libpaths=options.L)


    #buildModule('core/conformation', bindings_path, include_paths=options.I,
    #            libpaths=options.L, runtime_libpaths=options.L)


    updateList = [ (IncludeDictNew,'IncludeDict.new') ]

    if options.sort_IncludeDict:
        updateList.append( (IncludeDict,'IncludeDict') )

    for dc, fl in updateList:
        if Options.verbose: print 'Updating include dictionary file %s...' % fl

        def f_cmp(a, b):
            a_, b_ = a.split('/'), b.split('/')
            if cmp( a_[:-1], b_[:-1] ) : return cmp( a_[:-1], b_[:-1] )
            elif dc[a][1] == dc[b][1]: return cmp( a_[-1], b_[-1] )
            else: return cmp(dc[a][1], dc[b][1])


        K = dc.keys();  K.sort(cmp=f_cmp)
        f = file('python/bindings/' + fl + '.tmp', 'w');  f.write('#(True/False - include/exclude file/dir, priority, dependency)\n{\n')
        for k in K:
            #f.write( '%s : %s,\n' % (repr(k), repr({True:'+', False:'-'}[ dc[k] ]) ) )
            #if dc[k][2] == []: dc[k] = (dc[k][0], dc[k][1], None)
            f.write( '%s : %s,\n' % (repr(k), repr(dc[k]) ) )

        f.write('}\n');  f.close()
        os.rename('python/bindings/' + fl + '.tmp', 'python/bindings/' + fl)


    print "Done!"


def Sleep(time_, message, dict_={}):
    ''' Fancy sleep function '''
    len_ = 0
    for i in range(time_, 0, -1):
        #print "Waiting for a new revision:%s... Sleeping...%d     \r" % (sc.revision, i),
        msg = message % dict(dict_, time_left=i)
        print msg,
        len_ = max(len_, len(msg))
        sys.stdout.flush()
        time.sleep(1)

    print ' '*len_ + '\r',  # erazing sleep message


def SleepPrecise(time_, message, dict_={}):
    ''' Fancy sleep function '''
    if time_ > 0:
        len_ = 0
        msg = message % dict(dict_, time=time_)
        print msg,
        len_ = max(len_, len(msg))
        sys.stdout.flush()
        time.sleep(time_)

        print ' '*len_ + '\r',  # erazing sleep message


def execute(message, command_line, return_=False, untilSuccesses=False, print_output=True, verbose=True):
    if verbose:
        print message
        print command_line

    while True:
        #(res, output) = commands.getstatusoutput(commandline)

        po = subprocess.Popen(command_line+ ' 1>&2', bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #po = subprocess.Popen(command_line+ ' 1>&2', bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        f = po.stderr
        output = ''
        for line in f:
            #po.poll()
            if print_output: print line,
            output += line
            sys.stdout.flush()
        f.close()
        while po.returncode is None: po.wait()
        res = po.returncode
        #print '_____________________', res


        '''
        po = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderror = po.communicate()
        output = stdout + stderror
        res = po.returncode
        '''
        #print output

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if res:
        if print_output: print "\nEncounter error while executing: " + command_line
        if not return_: sys.exit(1)

    if return_ == 'output': return output
    else: return res


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
    if pid: Jobs.append( NT(pid=pid, tag=tag) )  # We are parent!
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
    add_option += (' -DDEBUG' if Options.debug else ' -DNDEBUG')

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
    else:
        add_loption = '-dynamiclib -Xlinker -headerpad_max_install_names'
        #if platform.release()[:2] == '13': add_loption += ' -stdlib=libstdc++'

    return add_loption




def buildModules__old(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
    ''' recursive build buinding for given dir name, and store them in dest.
    '''
    def visit(arg, dir_name, names):
        if dir_name.find('.svn') >= 0: return  # exclude all svn related namespaces

        # temp  for generating original exclude list
        #IncludeDict[dir_name] = not exclude.namespace(dir_name)

        #if exclude.namespace(dir_name): return

        if Options.build_all:
            if exclude.isBanned(dir_name):
                if Options.verbose: print 'Dir %s is banned! Skipping...' % dir_name
                return
        else:
            if dir_name in IncludeDict:
                if not IncludeDict[dir_name][0]:
                    print 'Skipping dir %s...' % dir_name
                    return
            else:
                print "Skipping new dir", dir_name
                IncludeDictNew[dir_name] = (False, 999, [])
                return



        #print "buildModules(...): '%s', " % dir_name
        #print "Directory: ", dir_name
        #dname = dest+'/' + os.path.dirname(dir_name)
        dname = dest+'/' + dir_name
        if not os.path.isdir(dname): os.makedirs(dname)
        '''
        if Options.parsing_jobs > 1:
            sys.stdout.flush()
            pid = mFork()
            if not pid:  # we are child process
                buildModule(dir_name, dest, include_paths, libpaths, runtime_libpaths, gccxml_path)
                sys.exit(0)

        else:
            IncludeDictNew.update( buildModule(dir_name, dest, include_paths, libpaths, runtime_libpaths, gccxml_path) )
            '''
        IncludeDictNew.update( buildModule(dir_name, dest, include_paths, libpaths, runtime_libpaths, gccxml_path) )

    os.path.walk(path, visit, None)


def buildModules(paths, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
    ''' recursive build buinding for given dir name, and store them in dest.
    '''
    #os.path.walk(path, visit, None)
    dir_list = []
    for path in paths:
        for dir_name, _, files in os.walk(path):
            if dir_name.find('.svn') >= 0: continue  # exclude all svn related namespaces

            if Options.build_all:
                if exclude.isBanned(dir_name):
                    if Options.verbose: print 'Dir %s is banned! Skipping...' % dir_name
                    continue
            else:
                if dir_name in IncludeDict:
                    if not IncludeDict[dir_name][0]:
                        print 'Skipping dir %s...' % dir_name
                        continue
                else:
                    print "Skipping new dir", dir_name
                    IncludeDictNew[dir_name] = (False, 999, [])
                    continue

            dir_list.append( (dir_name, files) )

    dir_list.sort(key=lambda x: -len(x[1]))  # sort dirs by number of files, most populated first. This should improve speed of multi-thread builds
    #for d, fs in dir_list: print len(fs), d

    if Options.one_lib_file and Options.build_all:
        mb = []
        for dir_name, _ in dir_list:
            #print "buildModules(...): '%s', " % dir_name
            dname = dest+'/' + dir_name
            if not os.path.isdir(dname): os.makedirs(dname)

            mb.append( ModuleBuilder(dir_name, dest, include_paths, libpaths, runtime_libpaths, gccxml_path) )
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

    else:
        for dir_name, _ in dir_list:
            print "buildModules(...): '%s', " % dir_name
            #print "Directory: ", dir_name
            #dname = dest+'/' + os.path.dirname(dir_name)
            dname = dest+'/' + dir_name
            if not os.path.isdir(dname): os.makedirs(dname)

            IncludeDictNew.update( buildModule(dir_name, dest, include_paths, libpaths, runtime_libpaths, gccxml_path) )





def prepareMiniLibs(mini_path, bindings_path):
    mode = 'pyrosetta_debug' if Options.debug else 'pyrosetta'

    if Platform == "macos" and PlatformBits=='32': execute("Building Rosetta libraries...", "cd %s && ./scons.py mode=%s arch=x86 arch_size=32 -j%s" % (mini_path, mode, Options.jobs) )
    elif Platform == "macos" and PlatformBits=='64': execute("Building mini libraries...", "cd %s && ./scons.py mode=%s -j%s" % (mini_path, mode, Options.jobs) )
    elif Platform == "cygwin": execute("Building mini libraries...", "cd %s && ./scons.py mode=%s bin -j%s" % (mini_path, mode, Options.jobs) )
    else: execute("Building mini libraries...", "cd %s && ./scons.py mode=%s -j%s" % (mini_path, mode, Options.jobs) )

    version_add_on = execute("Getting GCC version...", 'gcc -dumpversion', return_='output').strip()[0:3] + '/default/'

    # fix this for diferent platform
    if Platform == "linux": lib_path = 'build/src/'+mode+'/linux/' + platform.release()[:3] + '/' + PlatformBits +'/x86/gcc/'
    elif Platform == "cygwin": lib_path = 'build/src/'+mode+'/cygwin/1.7/32/x86/gcc/'
    else:
        if Platform == "macos" and PlatformBits=='32': lib_path = 'build/src/'+mode+'/macos/10.5/32/x86/gcc/'
        if Platform == "macos" and PlatformBits=='64':
            if platform.release()[:2] == '10': lib_path = 'build/src/'+mode+'/macos/10.6/64/x86/gcc/'
            elif platform.release()[:2] == '11': lib_path = 'build/src/'+mode+'/macos/10.7/64/x86/gcc/'
            elif platform.release()[:2] == '12': lib_path = 'build/src/'+mode+'/macos/10.8/64/x86/gcc/'
            elif platform.release()[:2] == '13': lib_path = 'build/src/'+mode+'/macos/10.9/64/x86/clang/';  version_add_on = '5.0/default/'
            else: print 'Unknown MacOS version:', platform.release()[:2];  sys.exit(1)  #lib_path = 'build/src/'+mode+'/macos/10.8/64/x86/gcc/'

    # now lets add version to lib_path...
    lib_path += version_add_on

        #if Platform == "macos" and PlatformBits=='64'  and  platform.release().startswith('11.'): lib_path = 'build/src/pyrosetta/macos/11/64/x86/gcc/'
        #else: lib_path = 'build/src/pyrosetta/macos/10.6/64/x86/gcc/'

    #lib_path += 'static/'
    obj_suffix = '.os'
    #if Platform == "linux"  and  PlatformBits == '64': obj_suffix = '.os'

    # Now the funny part - we rebuild all libs to produce just one lib file...
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
    if Platform == 'cygwin' : suffix = 'dll'
    if Platform == 'macos' : suffix = 'dylib'

    mini = mini_path + '/src/python/bindings/rosetta/libmini.' + suffix

    add_loption = getLinkerOptions()

    execute("Linking mini lib...",
            "cd %(mini_path)s && cd %(lib_path)s && gcc %(add_loption)s \
            %(objs)s -lz -lstdc++ -o %(mini)s" % dict(mini_path=mini_path, lib_path=lib_path, add_loption=add_loption, mini=mini, objs=objs, compiler=Options.compiler)
             )

    #if Platform == 'cygwin':
    #        execute("cp libs...", "cd %s && cp %s/*.dll %s" % (mini_path, lib_path, bindings_path) )
    #else:
    #        execute("cp libs...", "cd %s && cp %s/lib* %s" % (mini_path, lib_path, bindings_path) )

    if Platform == 'macos':
        #libs = ['libObjexxFCL.dylib', 'libnumeric.dylib', 'libprotocols.dylib', 'libdevel.dylib', 'libutility.dylib', 'libcore.dylib']
        libs = ['libmini.dylib']
        for l in libs:
            execute('Adjustin lib self path in %s' % l, 'install_name_tool -id rosetta/%s %s' % (l, bindings_path+'/'+l) )
            for k in libs:
                execute('Adjustin lib path in %s' % l, 'install_name_tool -change %s rosetta/%s %s' % (os.path.abspath(mini_path+'/'+lib_path+k), k, bindings_path+'/'+l) )

    shutil.copyfile(bindings_path+'/../src/__init__.py' , bindings_path+'/__init__.py')


def getAllRosettaSourceFiles():
    ''' return two lists: ([external] [rosetta]) source files '''
    obj_suffix = ''

    all_sources = []
    all_scons_files = [f for f in os.listdir('./../../')
        if f.endswith('.src.settings') and f not in ['apps.src.settings', 'pilot_apps.src.settings', 'devel.src.settings',]
           #and not f.startswith('protocols.')
        ]

    for scons_file in all_scons_files:
        all_sources.append([])

        f = file('./../../%s' % scons_file).read();  exec(f)
        for k in sources:
            for f in sources[k]:
                #all_sources.append( scons_file + '/' + k + '/' + f + obj_suffix)
                if not f.endswith('.cu'): all_sources[-1].append( k + '/' + f + '.cc')

    extra_objs = [  # additioanal source for external libs
        'dbio/cppdb/atomic_counter.cpp', "dbio/cppdb/conn_manager.cpp", "dbio/cppdb/driver_manager.cpp", "dbio/cppdb/frontend.cpp",
        "dbio/cppdb/backend.cpp", "dbio/cppdb/mutex.cpp", "dbio/cppdb/pool.cpp", "dbio/cppdb/shared_object.cpp", "dbio/cppdb/sqlite3_backend.cpp",
        "dbio/cppdb/utils.cpp", 'dbio/sqlite3/sqlite3.c', ]


    #all_sources.sort()  # for some reason sorting reduce lib size...
    #for i in all_sources: print i
    return extra_objs, all_sources


def get_CL_Options():
    if Options.debug:
        # Windows MSVC compiler common options (no optimization, but it works)
        return '/MD /GR /Gy /NDEBUG /DWIN32 /DBOOST_NO_MT /DPYROSETTA /DWIN_PYROSETTA /EHsc /nologo'  # /DNDEBUG
    else:
        # MSVC release options: /O2 /Oi /GL /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_UNICODE" /D "UNICODE" /FD /EHsc /MD /Gy /Yu"stdafx.h" /Fp"Release\SimpleConsoleApp.pch" /Fo"Release\\" /Fd"Release\vc90.pdb" /W3 /nologo /c /Zi /TP /errorReport:prompt
        # /Zi is 'generate omplete debugging information' - removing it
        # /W3 is just warning level
        # old one, 'the originale': return '/O2 /Oi /GL /DWIN32 /D "NDEBUG" /D "_CONSOLE" /FD /EHsc /MD /Gy /nologo /TP /DBOOST_NO_MT /DPYROSETTA /DWIN_PYROSETTA'
        # removing /GL
        # adding /bigobj
        return '/bigobj /O2 /Oi /DWIN32 /D "NDEBUG" /D "_CONSOLE" /FD /EHsc /MD /Gy /nologo /TP /DBOOST_NO_MT /DPYROSETTA /DWIN_PYROSETTA'



def BuildRosettaOnWindows(build_dir, bindings_path):
    ''' bypassing scones and build rosetta on windows native
    '''
    # copy svn_version file from pre-generated sources
    shutil.copy2(Options.use_pre_generated_sources + '/svn_version.cc', './../../core/svn_version.cc')

    external, sources = getAllRosettaSourceFiles()
    #if 'protocols/moves/PyMolMover.cc' in sources: sources.remove('protocols/moves/PyMolMover.cc')

    os.chdir( './../../' )

    #if Options.BuildMiniLibs:
    print 'Building Rosetta lib...'
    # Generate svn_version
    if (not os.path.isfile('core/svn_version.cc')): execute('Generate svn_version.cc...', 'cd .. && python svn_version.py')

    for s in external:
        try: os.makedirs( os.path.join( build_dir, os.path.split(s)[0]) )
        except OSError: pass

        obj = os.path.join(build_dir,s) + '.obj'
        s = '../external/' + s

        if (not os.path.isfile(obj))   or  os.path.getmtime(obj) < os.path.getmtime(s):
            execute('Compiling %s' % s, 'cl /MD /GR /Gy /D "WIN32" /D "_CONSOLE" /D "_UNICODE" /D "UNICODE" /DSQLITE_DISABLE_LFS \
                    /DSQLITE_OMIT_LOAD_EXTENSION /DSQLITE_THREADSAFE=0 /DCPPDB_EXPORTS /DCPPDB_LIBRARY_SUFFIX=\\".dylib\\" \
                    /DCPPDB_LIBRARY_PREFIX=\\"lib\\" /DCPPDB_DISABLE_SHARED_OBJECT_LOADING /DCPPDB_DISABLE_THREAD_SAFETY \
                    /DCPPDB_SOVERSION=\\"0\\" /DCPPDB_MAJOR=0 /DCPPDB_MINOR=3 /DCPPDB_PATCH=0 /DCPPDB_VERSION=\\"0.3.0\\"  \
                    /DCPPDB_WITH_SQLITE3 /DBOOST_NO_MT /DPYROSETTA /DWIN_PYROSETTA /DNDEBUG /c %s /I. /I../external/include /I../external/boost_1_46_1 \
                    /I../external/dbio /I../external/dbio/sqlite3 /Iplatform/windows/PyRosetta /Fo%s /EHsc' % (s, obj),)

        '''DNDEBUG -DCPPDB_EXPORTS -DCPPDB_LIBRARY_SUFFIX=\".dylib\" -DCPPDB_LIBRARY_PREFIX=\"lib\" -DCPPDB_DISABLE_SHARED_OBJECT_LOADING
           -DCPPDB_DISABLE_THREAD_SAFETY -DCPPDB_SOVERSION=\"0\" -DCPPDB_WITH_SQLITE3
        '''

    latest = None
    for s in sum(sources, []):
        try: os.makedirs( os.path.join( build_dir, os.path.split(s)[0]) )
        except OSError: pass

        obj = os.path.join(build_dir,s) + '.obj'

        if (not os.path.isfile(obj))   or  os.path.getmtime(obj) < os.path.getmtime(s):
            execute('Compiling %s' % s, 'cl %s /c %s /I. /I../external/include /I../external/boost_1_46_1 /I../external/dbio /Iplatform/windows/PyRosetta /Fo%s' % (get_CL_Options(), s, obj),
                    )  # return_=True)

        latest = max(latest, os.path.getmtime(obj) )


    sources = [external] + sources;  # After this moment there is no point of trating 'external' as special...
    #sources = [ sum(sources, []) ];  sources[0].sort()
    objs_all = ' '.join( [f + '.obj' for f in sum(sources, [])] )
    file(os.path.join(build_dir,'objs_all') , 'w').write( objs_all )


    # Now creating DLL
    dll = os.path.join(bindings_path, '..\\rosetta.dll')
    ''' There is no longer need for intermediate libs... but we will need it for PilotApp builds
    rosetta_lib_ = os.path.join(build_dir, 'rosetta_lib-%02d.lib')  # we purposly name 'real' lib somewhat different to avoid confusion with rosetta.lib file that belong to DLL
    rosetta_libs = [rosetta_lib_ % l for l in range(len(sources)) ]

    for i, rosetta_lib in enumerate(rosetta_libs):
        #objs = ' '.join( [ f + '.obj' for f in external] ) + ' ' + ' '.join( [f + '.obj' for f in sources[i]] )
        objs = ' '.join( [f + '.obj' for f in sources[i]] )

        file(os.path.join(build_dir,'objs-%02d' % i) , 'w').write( objs )
        #execute('Creating DLL %s...' % dll, 'cd %s && link /OPT:NOREF /dll @objs ..\\..\\external\\lib\\win_pyrosetta_z.lib /out:%s' % (build_dir, dll) )

        if (not os.path.isfile(rosetta_lib))   or  os.path.getmtime(rosetta_lib) < latest:
            # /INCREMENTAL:NO /LTCG
            execute('Creating lib %s...' % rosetta_lib, 'cd %s && lib @objs-%02d ..\\..\\external\\lib\\win_pyrosetta_z.lib Ws2_32.lib /out:%s' % (build_dir, i, rosetta_lib))
            latest = max(latest, os.path.getmtime(rosetta_lib) )
    '''
    # libcmt.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
    #execute('Creating DLL %s...' % s, 'cd ..\\build\\windows && cl /link /DLL @objs /OUT:%s' % dll )
    '''
    rosetta_lib = ' '.join( [ l for l in rosetta_libs] )
    pdb_test = 'apps/pilot/sergey/PDBTest.cc'
    pdb_test_exe = os.path.join(build_dir, 'PDBTest.exe')
    pdb_test_obj = os.path.join(build_dir, 'PDBTest.obj')
    if (not os.path.isfile(pdb_test_exe))   or  os.path.getmtime(pdb_test_exe) < os.path.getmtime(pdb_test)  or  os.path.getmtime(pdb_test_exe) < latest:
        #execute('Compiling test executable %s' % pdb_test, 'cl /c %s %s /I. /I../external/include /I../external/boost_1_46_1 /I../external/dbio /Iplatform/windows/PyRosetta %s /Fe%s /Fo%s' % (get_CL_Options(), pdb_test, rosetta_lib, pdb_test_exe, pdb_test_obj) )
        execute('Compiling test executable %s' % pdb_test, 'cl /c %s %s /I. /I../external/include /I../external/boost_1_46_1 /I../external/dbio /Iplatform/windows/PyRosetta /Fe%s /Fo%s' % (get_CL_Options(), pdb_test, pdb_test_exe, pdb_test_obj) )
        execute('Linking test executable %s' % pdb_test, 'link %s %s /out:%s' % (pdb_test_obj, rosetta_lib, pdb_test_exe) )

    '''

    print 'Building bindings...'

    #def visit(base_dir, dir_name, names):
    #    wn_buildOneNamespace(base_dir, dir_name, names, bindings_path, build_dir)
    #
    #os.path.walk(Options.use_pre_generated_sources, visit, Options.use_pre_generated_sources)

    objs = []
    symbols = []
    #latest = None  # keeping track of dates of local .def files
    for dir_name, _, files in os.walk(Options.use_pre_generated_sources):
        l = wn_buildOneNamespace(Options.use_pre_generated_sources, dir_name, files, bindings_path, build_dir, symbols)
        latest = max(l, latest)

        py_objs = os.path.join(build_dir,'py_objs' )
        f = file(py_objs, 'w');  f.write( ' '.join(objs) );  f.close()

                #print '____ Adding %s <-- %s' % (symbols[-1], l)

    symbols = list( set(symbols) )
    #file('._all_needed_symbols_', 'w').write(res)

    def_file = os.path.join(build_dir, 'rosetta_symbols.def')

    if (not os.path.isfile(dll))  or  os.path.getmtime(dll) < latest:
        print '\n\nWriting final export list... %s symbols...' % len(symbols)
        f = file(def_file, 'w'); f .write('LIBRARY rosetta\nEXPORTS\n  ' + '\n  '.join(symbols) + '\n' );  f.close()
        execute('Creating DLL %s...' % dll, 'cd %s && link /OPT:NOREF /dll @objs_all ..\\..\\external\\lib\\win_pyrosetta_z.lib Ws2_32.lib /DEF:%s /out:%s' % (build_dir, def_file, dll) )

    for dir_name, _, files in os.walk(Options.use_pre_generated_sources):
        wn_buildOneNamespace(Options.use_pre_generated_sources, dir_name, files, bindings_path, build_dir, link=True)

    print 'Done building PyRosetta bindings for Windows!'


    #for o in objs: print o

    #"c:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvars32.bat"

    # "public: virtual bool __thiscall ObjexxFCL::IndexRange::contains(class ObjexxFCL::IndexRange const &)const " (?contains@IndexRange@ObjexxFCL@@UBE_NABV12@@Z)

def wn_buildOneNamespace(base_dir, dir_name, files, bindings_path, build_dir, all_symbols=[], link=False):
    files = sorted( filter(lambda f: f.endswith('.cpp'), files) )
    sub_dir = dir_name[ len(base_dir)+1: ]

    obj_dir = os.path.join(bindings_path, sub_dir)
    if not os.path.isdir(obj_dir): os.makedirs(obj_dir)

    init_file_src = os.path.join(dir_name, '__init__.py')
    if os.path.isfile(init_file_src): shutil.copyfile(init_file_src, os.path.join(os.path.join(bindings_path, sub_dir), '__init__.py') )

    pyd = os.path.join( obj_dir, '__%s_all_at_once_.pyd' % dir_name.split('\\')[-1])
    symbols_file = os.path.join( obj_dir, 'symbols')  # list of symbols needed for this DLL, one per line

    latest = None

    if files:
        rosetta_lib = os.path.join(bindings_path, '..\\rosetta.lib')  # this is actually link to DLL, don't get confused it with rosettta_lib-%d.lib
        #rosetta_dll = os.path.join(bindings_path, '..\\rosetta.dll')
        #rosetta_lib = ' '.join( [ l for l in rosetta_libs] )

        #if (not os.path.isfile(pyd))   or  os.path.getmtime(pyd) < max( [os.path.getmtime( os.path.join(dir_name, f) ) for f in files] ):

        objs = []
        for f in files:
            source = os.path.join(dir_name, f)
            obj = os.path.join( obj_dir, f[:-3]+'obj')
            objs.append(obj)

            if (not os.path.isfile(obj))   or  os.path.getmtime(obj) < os.path.getmtime(source):
                execute('Compiling %s\\%s' % (dir_name, f), ('cl %s /c %s' % (get_CL_Options(), source)  #  /D__PYROSETTA_ONE_LIB__
                   + ' /I. /I../external/include /IC:/WPyRosetta/boost_1_47_0 /I../external/dbio /Iplatform/windows/PyRosetta'
                   + ' /Ic:\Python27\include /DWIN_PYROSETTA_PASS_2 /DBOOST_PYTHON_MAX_ARITY=25'
                   + ' /Fo%s ' % obj ) )
                #  /Iplatform/windows/32/msvc
                #  /I../external/boost_1_46_1
                #  /IBOOST_MSVC    /link rosetta_lib

                # c:\\mingw\\bin\\
                """execute('Compiling %s' % (dir_name+f), 'gcc -DPYROSETTA -c %s -I. \
                        -I../external/include -IC:/WPyRosetta/boost_1_47_0 -I../external/dbio -Iplatform/windows/PyRosetta \
    -Ic:\Python27\include -c -pipe -O3 -ffast-math -funroll-loops -finline-functions -DBOOST_PYTHON_MAX_ARITY=25 \
        -o %s' % (source, obj) )"""

            latest = max(latest, os.path.getmtime(obj) )

        #
        if (not os.path.isfile(symbols_file))   or  os.path.getmtime(symbols_file) < latest:

            dummy = os.path.join( obj_dir, '_dummy_') #_rosetta_.pyd' )

            #if (not os.path.isfile(pyd))   or  os.path.getmtime(pyd) < latest:
            res = execute('Getting list of missed symbols... Creating DLL %s...' % dummy,
                            'link %s ..\\external\\lib\\win_pyrosetta_z.lib /INCREMENTAL:NO /dll /libpath:c:/Python27/libs /libpath:c:/WPyRosetta/boost_1_47_0/stage/lib /out:%s' % (' '.join(objs), dummy), return_='output', print_output=False)
            symbols = []
            for l in res.split('\n'):
                if   l.find('error LNK2001: unresolved external symbol __DllMainCRTStartup') >=0 : continue  # error LNK2001: unresolved external symbol __DllMainCRTStartup@12
                elif l.find('error LNK2019: unresolved external symbol') >=0 : symbols.append( l.partition('" (')[2].partition(') referenced in function')[0] )
                elif l.find('error LNK2001: unresolved external symbol') >=0 : symbols.append( l.partition('" (')[2][:-2]) # + ' DATA')
            f = file(symbols_file, 'w');  f.write( '\n'.join(symbols) );  f.close()
            print '\nAdding %s symbols... Tottal now is:%s\n' % (len(symbols), len(set(all_symbols+symbols)))

        all_symbols.extend( file(symbols_file).read().split('\n') )

        if link:
            if (not os.path.isfile(pyd))   or  os.path.getmtime(pyd) < latest  or  os.path.getmtime(pyd) < os.path.getmtime(rosetta_lib):
                execute('Creating DLL %s...' % pyd,
                        'link  /INCREMENTAL:NO /dll /libpath:c:/Python27/libs /libpath:c:/WPyRosetta/boost_1_47_0/stage/lib  %s %s ..\\external\\lib\\win_pyrosetta_z.lib /out:%s' % (' '.join(objs), rosetta_lib, pyd) )
                #map(os.remove, objs)

    return latest



#-c -pipe -O3 -ffast-math -funroll-loops -finline-functions -fPIC -DBOOST_PYTHON_MAX_ARITY=25 -I../external/include  -I../external/dbio
#-I/Users/sergey/work/trunk/PyRosetta.develop.Python-2.7/PyRosetta.Develop.64/include -I/Users/sergey/work/trunk/PyRosetta.develop.Python-2.7/PyRosetta.Develop.64/include/boost
#-I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I../src/platform/linux -I../src


        #cl /MD p_qwe.cpp /EHsc /I. /Ic:\T\boost_1_47_0_ /Ic:\Python27\include /c /GR /Gy /W3 /GS-

    #global __global;  __global += 1
    #if __global>3: sys.exit(1)



_TestInludes_ = ''

def buildModule(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
    ''' Build one namespace and return dict of newly found heades.
    '''
    '''
    testName = 'python/bindings/TestIncludes.py'
    global _TestInludes_
    _TestInludes_ += 'import rosetta.%s\n' % path.replace('/', '.')
    f = file(testName, 'w');  f.write(_TestInludes_);  f.close()
    '''
    gc.collect()

    if Options.py_plus_plus:
        #if path == 'core':
        #    buildModule_All(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path)
        #else: buildModule_One(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path)
        return buildModule_One(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path)

    else:
        #print 'Using Clang...'
        return buildModule_UsingCppParser(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path)



class ModuleBuilder:
    def __init__(self, path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
        ''' Non recursive build buinding for given dir name, and store them in dest.
            path - relative path to namespace
            dest - path to root file destination, actual dest will be dest + path

            return dict of a newly found headers.

            This is a class because we want to generate path/name only once etc.
        '''
        global Options
        if Options.verbose: print 'CppXML route: ModuleBuilder.init...', path, dest

        self.path = path
        self.dest = dest


        # Creating list of headers
        self.headers = [os.path.join(path, d)
                            for d in os.listdir(path)
                                if os.path.isfile( os.path.join(path, d) )
                                    and d.endswith('.hh')
                                    and not d.endswith('.fwd.hh')
                                    ]

        self.headers.sort()

        for h in self.headers[:]:
            if exclude.isBanned(h):
                if Options.verbose: print "Banning header:", h
                self.headers.remove(h)

        if Options.verbose: print_(self.headers, color='black', bright=True)

        self.fname_base = dest + '/' + path

        if not os.path.isdir(self.fname_base): os.makedirs(self.fname_base)

        #print 'Creating __init__.py file...'
        f = file( dest + '/' + path + '/__init__.py', 'w');  f.close()
        if not self.headers:  return   # if source files is empty then __init__.py should be empty too

        def finalize_init(current_fname):
                #print 'Finalizing Creating __init__.py file...'
                f = file( dest + '/' + path + '/__init__.py', 'a');
                f.write('from %s import *\n' % os.path.basename(current_fname)[:-3]);
                f.close()

        self.finalize_init = finalize_init

        self.include_paths = ' -I'.join( [''] + include_paths + ['../src/platform/linux' , '../src'] )

        self.libpaths = ' -L'.join( ['', dest] + libpaths )
        self.runtime_libpaths = ' -Xlinker -rpath '.join( [''] + runtime_libpaths + ['rosetta'] )

        self.cpp_defines = '-DPYROSETTA -DPYROSETTA_DISABLE_LCAST_COMPILE_TIME_CHECK -DBOOST_NO_MT'
        if Options.cross_compile: self.cpp_defines += ' -I../src/platform/windows/PyRosetta'
        elif Platform == "macos": self.cpp_defines += ' -I../src/platform/macos'
        else: self.cpp_defines += ' -I../src/platform/linux'

        #self.gccxml_options = '--gccxml-compiler llvm-g++-4.2 -march=nocona' if Platform == "macos" else ''
        self.gccxml_options = ''
        if Options.gccxml_compiler: self.gccxml_options += '--gccxml-compiler ' + Options.gccxml_compiler
        elif Platform == 'macos':
            if platform.release()[:2] == '13': self.gccxml_options += '--gccxml-compiler gcc'  #  --gccxml-cxxflags "-stdlib=libstdc++"
            else: self.gccxml_options += '--gccxml-compiler llvm-g++-4.2'
        if Platform == 'macos': self.gccxml_options += ' -march=nocona'

        self.cc_files = []
        self.add_option  = getCompilerOptions()
        self.add_loption = getLinkerOptions()

        self.by_hand_beginning_file = dest+ '/../src/' + path + '/_%s__by_hand_beginning.cc' % path.split('/')[-1]
        self.by_hand_ending_file = dest+ '/../src/' + path + '/_%s__by_hand_ending.cc' % path.split('/')[-1]

        self.by_hand_beginning = file(self.by_hand_beginning_file).read() if os.path.isfile(self.by_hand_beginning_file) else ''
        self.by_hand_ending = file(self.by_hand_ending_file).read() if os.path.isfile(self.by_hand_ending_file) else ''

        self.all_at_once_base = '__' + path.split('/')[-1] + '_all_at_once_'
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
        namespace = os.path.basename(self.path)
        py_init_file = self.dest + '/../src/' + self.path + '/__init__.py'
        if os.path.isfile(py_init_file): t = file(py_init_file).read()
        else: t = ''
        f = file( self.dest + '/' + self.path + '/__init__.py', 'w');  f.write(t+'from %s import *\n' % self.all_at_once_base);  f.close()

        if xml_recompile or (not Options.update):
            start_time = time.time()

            if os.path.isfile(self.all_at_once_lib): os.remove(self.all_at_once_lib)

            def generate():
                if execute('Generating XML representation...', 'gccxml %s %s -fxml=%s %s -I. -I../external/include -I../external/boost_1_46_1  -I../external/dbio -DBOOST_NO_INITIALIZER_LISTS ' % (self.gccxml_options, self.all_at_once_source_cpp, self.all_at_once_xml, self.cpp_defines), Options.continue_ ): return

                namespaces_to_wrap = ['::'+self.path.replace('/', '::')+'::']
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

                    exclude.finalize2(all_at_once_N_cpp, self.dest, self.path, module_name=self.all_at_once_base, add_by_hand = False, includes=includes)
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
                SleepPrecise( (time.time() - start_time) * Options.sleep, 'Sleeping %(time)s seconds so CPU could cool-off...\r')


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
                comiler_cmd = "%(compiler)s %(fname)s -o %(obj_name)s -c %(add_option)s %(cpp_defines)s -I../external/include  -I../external/dbio %(include_paths)s "
                comiler_dict = dict(add_option=self.add_option, fname=all_at_once_N_cpp, obj_name=all_at_once_N_obj, include_paths=self.include_paths, compiler=Options.compiler, cpp_defines=self.cpp_defines)

                failed = False

                if not Options.cross_compile:
                    def compile_():
                        if execute("Compiling...", comiler_cmd % comiler_dict, return_=True):
                            if Options.compiler != 'clang': failed = True
                            elif execute("Compiling...", comiler_cmd % dict(comiler_dict, compiler='gcc'), return_=True): failed = True

                    if Options.jobs > 1:
                        pid = mFork(tag=self.path)
                        if not pid:  # we are child process
                            compile_()
                            sys.exit(0)

                    else:
                        compile_()
                        SleepPrecise( (time.time() - start_time) * Options.sleep, 'Sleeping %(time)s seconds so CPU could cool-off...\r')


                if Options.jobs == 1:
                    if Options.continue_ and failed: return new_headers


    def linkBindings(self):
        ''' Build early generated bindings.
        '''
        if not self.headers  or  Options.cross_compile: return
        source_list = json.load( file(self.all_at_once_json) )

        relink = False

        if Options.update:
            for (all_at_once_N_cpp, all_at_once_N_obj) in source_list:
                if not os.path.isfile(self.all_at_once_lib)  or  os.path.getmtime( all_at_once_N_obj) > os.path.getmtime(self.all_at_once_lib): relink = True; break
        else: relink = True


        if relink:
            if not Options.cross_compile:  # -fPIC -ffloat-store -ffor-scope
                start_time = time.time()

                objs_list = map(lambda x:x[1], source_list)
                linker_cmd = "cd %(dest)s/../ && %(compiler)s %(obj)s %(add_option)s -lmini -lstdc++ -lz -l%(python_lib)s \
                              -l%(boost_lib)s %(libpaths)s %(runtime_libpaths)s -o %(dst)s"
                linker_dict = dict(add_option=self.add_loption, obj=' '.join(objs_list), dst=self.all_at_once_lib, libpaths=self.libpaths, runtime_libpaths=self.runtime_libpaths, dest=self.dest, boost_lib=Options.boost_lib,
                        python_lib=Options.python_lib, compiler=Options.compiler)
                def linking():
                    if execute("Linking...", linker_cmd % linker_dict, return_= (True if Options.compiler != 'gcc' or Options.continue_ else False) ):
                        if Options.compiler != 'gcc':
                            execute("Linking...", linker_cmd % dict(linker_dict, compiler='gcc'), return_= Options.continue_ )

                if Options.jobs > 1:

                    pid = mFork(tag=self.path+'+linking', overhead=1)  # we most likely can start extra linking process, beceause it depend on compilation to  finish. There is no point of waiting for it...
                    if not pid:  # we are child process
                        #mWait(tag=self.path)  # wait for all compilation jobs to finish...
                        linking()
                        sys.exit(0)
                else:
                    linking()
                    SleepPrecise( (time.time() - start_time) * Options.sleep, 'Sleeping %(time)s seconds so CPU could cool-off...\r')


            else: execute("Toching %s file..." % self.all_at_once_lib, 'cd %(dest)s/../ && touch %(dst)s' % dict(dest=self.dest, dst=self.all_at_once_lib) )

        #print 'Done!'



def buildModule_UsingCppParser(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
    ''' Non recursive build buinding for given dir name, and store them in dest.
        path - relative path to namespace
        dest - path to root file destination, actual dest will be dest + path

        return dict of a newly found headers.
    '''
    print 'CppXML route: Processing', path, dest

    global Options

    # Creating list of headers
    if os.path.isfile(path): headers = [path]
    else:  headers = [os.path.join(path, d)
                        for d in os.listdir(path)
                            if os.path.isfile( os.path.join(path, d) )
                                and d.endswith('.hh')
                                and not d.endswith('.fwd.hh')
                                ]

    #headers.sort()
    # now, because of abstract class issue we have to sort header list first...
    #headers.sort(key = lambda x: IncludeDict.get(x, (None, 999, None) )[1])
    headers.sort()

    # tmp  for generating original exclude list
    #for i in headers: IncludeDict[i] = True

    new_headers = {}
    for h in headers[:]:

        if Options.build_all:
            if exclude.isBanned(h):
                print "Banning header:", h
                headers.remove(h)

            continue  # do not exclude anything...

        if h in IncludeDict:
            #print 'IncludeDict[%s] --> %s' % (h, IncludeDict[h])
            if not IncludeDict[h][0]:
                print "Excluding header:", h
                headers.remove(h)

        else:
            print "Excluding new header:", h
            headers.remove(h)
            new_headers[h] = (False, 999, [])

    # Temporary injecting Mover in to protocols level
    #if path == 'protocols': headers.insert(0, 'protocols/moves/Mover.hh')

    print headers

    fname_base = dest + '/' + path

    if not os.path.isdir(fname_base): os.makedirs(fname_base)

    print 'Creating __init__.py file...'
    f = file( dest + '/' + path + '/__init__.py', 'w');  f.close()
    if not headers:  return new_headers  # if source files is empty then __init__.py should be empty too
    def finalize_init(current_fname):
            #print 'Finalizing Creating __init__.py file...'
            f = file( dest + '/' + path + '/__init__.py', 'a');
            f.write('from %s import *\n' % os.path.basename(current_fname)[:-3]);
            f.close()


    include_paths = ' -I'.join( [''] + include_paths + ['../src/platform/linux' , '../src'] )

    nsplit = 1  # 1 files per iteration

    libpaths = ' -L'.join( ['', dest] + libpaths )
    runtime_libpaths = ' -Xlinker -rpath '.join( [''] + runtime_libpaths + ['rosetta'] )
    #if Platform == 'cygwin': runtime_libpaths = ' '

    cpp_defines = '-DPYROSETTA -DPYROSETTA_DISABLE_LCAST_COMPILE_TIME_CHECK -DBOOST_NO_MT'
    if Options.cross_compile: cpp_defines += ' -I../src/platform/windows/PyRosetta'
    elif Platform == "macos": cpp_defines += ' -I../src/platform/macos'
    else: cpp_defines += ' -I../src/platform/linux'

    gccxml_options = '--gccxml-compiler llvm-g++-4.2 -march=nocona' if Platform == "macos" else ''

    cc_files = []
    add_option  = getCompilerOptions()
    #if Options.compiler != 'gcc': add_option += ' -Wno-local-type-template-args'
    add_loption = getLinkerOptions()

    #by_hand_code = dest+ '/../src/' + path + '/_%s__by_hand.cc' % path.split('/')[-1]
    by_hand_beginning_file = dest+ '/../src/' + path + '/_%s__by_hand_beginning.cc' % path.split('/')[-1]
    by_hand_ending_file = dest+ '/../src/' + path + '/_%s__by_hand_ending.cc' % path.split('/')[-1]

    by_hand_beginning = file(by_hand_beginning_file).read() if os.path.isfile(by_hand_beginning_file) else ''
    by_hand_ending = file(by_hand_ending_file).read() if os.path.isfile(by_hand_ending_file) else ''

    #print '_______________ by_hand_ending=', by_hand_ending

    if Options.one_lib_file:
        all_at_once_base = '__' + path.split('/')[-1] + '_all_at_once_'
        all_at_once_source_cpp = fname_base + '/' + all_at_once_base + '.source.cc'
        all_at_once_cpp = fname_base + '/' + all_at_once_base + '.'
        all_at_once_obj = fname_base + '/' + all_at_once_base + '.'
        all_at_once_xml = fname_base + '/' + all_at_once_base + '.xml'
        all_at_once_lib = fname_base + '/' + all_at_once_base + '.so'
        if Platform == 'cygwin': all_at_once_lib = all_at_once_lib = fname_base + '/' + all_at_once_base + '.dll'
        all_at_once_relative_files = []
        xml_recompile = False

    for fl in headers:
        #print '\033[32m\033[1m%s\033[0m' % fl
        print_(fl, color='green', bright=True)

        #print 'Binding:', files
        hbase = fl.split('/')[-1][:-3]
        hbase = hbase.replace('.', '_')
        #print 'hbase = ', hbase
        #if hbase == 'init': hbase='tint'  # for some reason Boost don't like 'init' name ? by hand?

        fname =     fname_base + '/' + '_' + hbase + '.cc'
        inc_name =  fname_base + '/' + '_' + hbase + '.hh'
        obj_name =  fname_base + '/' + '_' + hbase + '.o'
        xml_name =  fname_base + '/' + '_' + hbase + '.xml'
        cc_for_xml_name =  fname_base + '/' + '_' + hbase + '.xml.cpp'
        dst_name =  fname_base + '/' + '_' + hbase + '.so'
        if Platform == 'cygwin' : dst_name =  fname_base + '/' + '_' + hbase + '.dll'
        decl_name = fname_base + '/' + '_' + hbase + '.exposed_decl.pypp.txt'

        cc_files.append(fname)

        if Options.update:
            try:
                if fl == headers[0]:  # for first header we additionaly check if 'by_hand' code is up to date
                    #print '__checking for:', by_hand_code
                    if Options.one_lib_file:
                        #if os.path.isfile(by_hand_code) and os.path.getmtime(by_hand_code) > os.path.getmtime(all_at_once_lib):
                        #    raise os.error

                        if os.path.isfile(by_hand_beginning_file) and os.path.getmtime(by_hand_beginning_file) > os.path.getmtime(all_at_once_lib): raise os.error
                        if os.path.isfile(by_hand_ending_file) and os.path.getmtime(by_hand_ending_file) > os.path.getmtime(all_at_once_lib): raise os.error

                    else:
                        #if os.path.isfile(by_hand_code) and os.path.getmtime(by_hand_code) > os.path.getmtime(dst_name): raise os.error

                        if os.path.isfile(by_hand_beginning_file) and os.path.getmtime(by_hand_beginning_file) > os.path.getmtime(dst_name): raise os.error
                        if os.path.isfile(by_hand_ending_file) and os.path.getmtime(by_hand_ending_file) > os.path.getmtime(dst_name): raise os.error

                if Options.one_lib_file:
                    if os.path.getmtime(fl) > os.path.getmtime(all_at_once_lib):
                        #print 'xml_recompile = True'

                        xml_recompile = True
                    else:
                        print 'File: %s is up to date - skipping' % fl
                        #finalize_init(fname)
                        #continue
                else:
                    if os.path.getmtime(fl) < os.path.getmtime(dst_name):
                        print 'File: %s is up to date - skipping' % fl
                        finalize_init(fname)
                        continue
            except os.error: xml_recompile = True

        source_fwd_hh = fl.replace('.hh', '.fwd.hh')
        source_hh = fl
        source_cc = fl.replace('.hh', '.cc')

        if Options.one_lib_file:  # just collecting file names...
            all_at_once_relative_files.extend( [source_fwd_hh, source_hh, source_cc] )
            continue


        source_to_parse = fl.replace('.hh', '.cc')
        if not os.path.isfile( source_cc ):  # there is no cc file... ^_^ - generating it on the fly... (clang dont like .hh etc)
            _h = file(cc_for_xml_name, 'w');  _h.write('#include <%s>\n' % fl);  _h.close()
            source_cc = cc_for_xml_name
        # -cc1 -ast-print-xml
        #execute('pwd...', 'pwd' ) # -I../external/include  -I.
        #execute('Generating XML representation...', 'clang++ -cc1 -ast-print-xml %s  ' % (cc_original,) )

        # Clang++ version
        #execute('Generating XML representation...', 'clang++ -cc1 -ast-print-xml %s -o %s -I. -I../src/platform/linux  -I../external/include -I../external/boost_1_46_1 ' % (source_cc, xml_name) )

        # we need  -DBOOST_NO_INITIALIZER_LISTS or gccxml choke on protocols/genetic_algorithm/GeneticAlgorithm.hh
        # GCCXML version
        # -DPYROSETTA
        if execute('Generating XML representation...', 'gccxml %s %s -fxml=%s %s -I. -I../external/include -I../external/boost_1_46_1 -I../external/dbio -DBOOST_NO_INITIALIZER_LISTS ' % (gccxml_options, source_hh, xml_name, cpp_defines), Options.continue_): continue

        namespaces_to_wrap = ['::'+path.replace('/', '::')+'::']
        # Temporary injecting Mover in to protocols level
        #if path == 'protocols': namespaces_to_wrap.append('::protocols::moves::')

        code = tools.CppParser.parseAndWrapModule('_'+hbase, namespaces_to_wrap,  xml_name, [source_fwd_hh, source_hh, source_cc],
                                                  by_hand_beginning = by_hand_beginning if fl==headers[0] else '',
                                                  by_hand_ending= by_hand_ending if fl==headers[-1] else '')
        if len(code) != 1:
            print 'Whats going on???'
            sys.exit(1)
        else: code = code[0]

        #code = '#include <%s>\n' % fl + code

        f = file(fname, 'w');  f.write(code);  f.close()

        #mb.build_code_creator( module_name = '_'+hbase, doc_extractor=doxygen.doxygen_doc_extractor() )
        #mb.code_creator.user_defined_directories.append( os.path.abspath('.') )  # make include relative
        #mb.write_module( os.path.join( os.path.abspath('.'), fname ) )

        #exclude.finalize(fname, dest, path, None, module_name='_'+hbase, add_by_hand = (fl==headers[0]), files=[fl], add_includes=True)
        exclude.finalize(fname, dest, path, None, module_name='_'+hbase, add_by_hand = False, files=[fl], add_includes=True)
        finalize_init(fname)

        if not Options.one_lib_file:
            if execute("Compiling...", # -fPIC
                "%(compiler)s %(fname)s -o %(obj_name)s -c \
                 %(add_option)s %(cpp_defines)s -I../external/include -I../external/dbio \
                 %(include_paths)s " % dict(add_option=add_option, fname=fname, obj_name=obj_name, include_paths=include_paths, compiler=Options.compiler, cpp_defines=cpp_defines),
                 Options.continue_):
                pass
                '''
                print 'Compiliation failed... Creating empty C++ file to test includes...'

                fname =    fname_base + '/' + '_0_.cc'
                obj_name = fname_base + '/' + '_0_.o'
                f = file(fname, 'w');  f.write('#include <%s>\n' % fl);  f.close()
                execute("Compiling empty C++...",
                    "gcc %(fname)s -o %(obj_name)s -c \
                     %(add_option)s -I../external/include \
                     %(include_paths)s " % dict(add_option=add_option, fname=fname, obj_name=obj_name, include_paths=include_paths))
                '''
                #sys.exit(1)

            execute("Linking...", # -fPIC -ffloat-store -ffor-scope
                "cd %(dest)s/../ && %(compiler)s %(obj)s %(add_option)s  \
                -lmini -lstdc++ -lz \
                 -l%(python_lib)s \
                 -l%(boost_lib)s \
                 %(libpaths)s %(runtime_libpaths)s -o %(dst)s" % dict(add_option=add_loption, obj=obj_name,
                    dst=dst_name, libpaths=libpaths, runtime_libpaths=runtime_libpaths, dest=dest, boost_lib=Options.boost_lib,
                    python_lib=Options.python_lib, compiler=Options.compiler
                    ),
                 Options.continue_)

    if Options.one_lib_file:
        f = file(all_at_once_source_cpp, 'w');
        for fl in headers: f.write('#include <%s>\n' % fl);
        f.close()

        #hbase = all_at_once_cpp.split('/')[-1][:-3]
        #hbase = hbase.replace('.', '_')

        #print 'Finalizing Creating __init__.py file...'
        namespace = os.path.basename(path)
        py_init_file = dest + '/../src/' + path + '/__init__.py'
        if os.path.isfile(py_init_file): t = file(py_init_file).read()
        else: t = ''
        f = file( dest + '/' + path + '/__init__.py', 'w');  f.write(t+'from %s import *\n' % all_at_once_base);  f.close()

        if xml_recompile or (not Options.update):
            if os.path.isfile(all_at_once_lib): os.remove(all_at_once_lib)

            if execute('Generating XML representation...', 'gccxml %s %s -fxml=%s %s -I. -I../external/include -I../external/boost_1_46_1  -I../external/dbio -DBOOST_NO_INITIALIZER_LISTS ' % (gccxml_options, all_at_once_source_cpp, all_at_once_xml, cpp_defines), Options.continue_ ):
                return new_headers

            namespaces_to_wrap = ['::'+path.replace('/', '::')+'::']
            # Temporary injecting Mover in to protocols level
            #if path == 'protocols': namespaces_to_wrap.append('::protocols::moves::')

            code = tools.CppParser.parseAndWrapModule(all_at_once_base, namespaces_to_wrap,  all_at_once_xml, all_at_once_relative_files, max_funcion_size=Options.max_function_size,
                                                      by_hand_beginning=by_hand_beginning, by_hand_ending=by_hand_ending)

            objs_list = []

            print_('Getting include list...', color='black', bright=True)
            includes = exclude.getIncludes(headers)

            print_('Finalizing[%s]' % len(code), color='black', bright=True, endline=False);  sys.stdout.flush()
            for i in range( len(code) ):
                all_at_once_N_cpp = all_at_once_cpp+'%s.cpp' % i
                all_at_once_N_obj = all_at_once_obj+'%s.o' % i
                if os.path.isfile(all_at_once_N_obj): os.remove(all_at_once_N_obj)

                for fl in headers: code[i] = '#include <%s>\n' % fl + code[i]

                f = file(all_at_once_N_cpp, 'w');  f.write(code[i]);  f.close()

                exclude.finalize2(all_at_once_N_cpp, dest, path, module_name=all_at_once_base, add_by_hand = False, includes=includes)
                print_('.', color='black', bright=True, endline=False); sys.stdout.flush()
            print_(' Done!', color='black', bright=True);

            for i in range( len(code) ):
                all_at_once_N_cpp = all_at_once_cpp+'%s.cpp' % i
                all_at_once_N_obj = all_at_once_obj+'%s.o' % i

                # -fPIC
                comiler_cmd = "%(compiler)s %(fname)s -o %(obj_name)s -c %(add_option)s %(cpp_defines)s -I../external/include  -I../external/dbio %(include_paths)s "
                comiler_dict = dict(add_option=add_option, fname=all_at_once_N_cpp, obj_name=all_at_once_N_obj, include_paths=include_paths, compiler=Options.compiler, cpp_defines=cpp_defines)

                failed = False

                if not Options.cross_compile:
                    def compile_():
                        if execute("Compiling...", comiler_cmd % comiler_dict, return_=True):
                            if Options.compiler != 'clang': failed = True
                            elif execute("Compiling...", comiler_cmd % dict(comiler_dict, compiler='gcc'), return_=True): failed = True

                    if Options.jobs > 1:
                        pid = mFork(tag=path)
                        if not pid:  # we are child process
                            compile_()
                            sys.exit(0)

                    else:
                        compile_()

                if Options.jobs == 1:
                    if Options.continue_ and failed: return new_headers

                    objs_list.append(all_at_once_N_obj)

            if not Options.cross_compile:  # -fPIC -ffloat-store -ffor-scope

                linker_cmd = "cd %(dest)s/../ && %(compiler)s %(obj)s %(add_option)s -lmini -lstdc++ -lz -l%(python_lib)s \
                              -l%(boost_lib)s %(libpaths)s %(runtime_libpaths)s -o %(dst)s"
                linker_dict = dict(add_option=add_loption, obj=' '.join(objs_list), dst=all_at_once_lib, libpaths=libpaths, runtime_libpaths=runtime_libpaths, dest=dest, boost_lib=Options.boost_lib,
                        python_lib=Options.python_lib, compiler=Options.compiler)
                def linking():
                    if execute("Linking...", linker_cmd % linker_dict, return_= (True if Options.compiler != 'gcc' or Options.continue_ else False) ):
                        if Options.compiler != 'gcc':
                            execute("Linking...", linker_cmd % dict(linker_dict, compiler='gcc'), return_= Options.continue_ )

                if Options.jobs > 1:

                    pid = mFork(tag=path+'+linking', overhead=1)  # we most likely can start extra linking process, beceause it depend on compilation to  finish. There is no point of waiting for it...
                    if not pid:  # we are child process
                        mWait(tag=path)  # wait for all compilation jobs to finish...
                        linking()
                        sys.exit(0)
                else:
                    linking()


            else: execute("Toching %s file..." % all_at_once_lib, 'cd %(dest)s/../ && touch %(dst)s' % dict(dest=dest, dst=all_at_once_lib) )




    print 'Done!', new_headers
    return new_headers





def buildModule_One(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
    ''' Non recursive build buinding for given dir name, and store them in dest.
        path - relative path to namespace
        dest - path to root file destination, actual dest will be dest + path

        return dict of a newly found headers.
    '''

    #if ParallelBuild :  # When building in parallel we have to set path for each thread as new...
    #    os.chdir( './../../' )


    print 'Processing', path

    # Creating list of headers
    headers = [os.path.join(path, d)
                for d in os.listdir(path)
                    if os.path.isfile( os.path.join(path, d) )
                        and d.endswith('.hh')
                        and not d.endswith('.fwd.hh')
                        ]
    #headers.sort()
    # now, because of abstract class issue we have to sort header list first...
    headers.sort(key = lambda x: IncludeDict.get(x, (None, 999, None) )[1])

    # tmp  for generating original exclude list
    #for i in headers: IncludeDict[i] = True

    new_headers = {}
    for h in headers[:]:
        if h in IncludeDict:
            #print 'IncludeDict[%s] --> %s' % (h, IncludeDict[h])
            if not IncludeDict[h][0]:
                print "Excluding header:", h
                headers.remove(h)
            #else:
            #    print 'Including %s' % h
        else:
            print "Excluding new header:", h
            headers.remove(h)
            new_headers[h] = (False, 999, [])


    print headers

    fname_base = dest + '/' + path

    if not os.path.isdir(fname_base): os.makedirs(fname_base)

    print 'Creating __init__.py file...'
    f = file( dest + '/' + path + '/__init__.py', 'w');  f.close()
    if not headers:  return new_headers  # if source files is empty then __init__.py should be empty too
    def finalize_init(current_fname):
            #print 'Finalizing Creating __init__.py file...'
            f = file( dest + '/' + path + '/__init__.py', 'a');
            f.write('from %s import *\n' % os.path.basename(current_fname)[:-3]);
            f.close()


    include_paths = ' -I'.join( [''] + include_paths + ['../src/platform/linux' , '../src'] )

    nsplit = 1  # 1 files per iteration

    libpaths = ' -L'.join( ['', dest] + libpaths )
    runtime_libpaths = ' -Xlinker -rpath '.join( [''] + runtime_libpaths + ['rosetta'] )
    #if Platform == 'cygwin': runtime_libpaths = ' '

    global Options

    cc_files = []
    add_option  = getCompilerOptions()
    add_loption = getLinkerOptions()

    #last_header_list = []
    #last_header_list_ = []
    for fl in headers:
        #last_header_list.extend(last_header_list_)
        #last_header_list_ = [fl]

        #fcount += 1
        #files = headers[:nsplit]
        #headers = headers[nsplit:]
        #print 'Binding:', files
        hbase = fl.split('/')[-1][:-3]
        hbase = hbase.replace('.', '_')
        #print 'hbase = ', hbase
        #if hbase == 'init': hbase='tint'  # for some reason Boost don't like 'init' name ? by hand?

        fname =     fname_base + '/' + '_' + hbase + '.cc'
        inc_name =  fname_base + '/' + '_' + hbase + '.hh'
        obj_name =  fname_base + '/' + '_' + hbase + '.o'
        dst_name =  fname_base + '/' + '_' + hbase + '.so'
        if Platform == 'cygwin' : dst_name =  fname_base + '/' + '_' + hbase + '.dll'
        decl_name = fname_base + '/' + '_' + hbase + '.exposed_decl.pypp.txt'

        cc_files.append(fname)

        if Options.update:
            try:
                if Options.one_lib_file:
                    if os.path.getmtime(fl) < os.path.getmtime(fname):
                        print 'File: %s is up to date - skipping' % fl
                        finalize_init(fname)
                        continue

                else:
                    if os.path.getmtime(fl) < os.path.getmtime(dst_name):
                        print 'File: %s is up to date - skipping' % fl
                        finalize_init(fname)
                        continue
            except os.error: pass


        mb = module_builder.module_builder_t(files= [fl]
                                         , include_paths = ['../src/platform/linux', '../external/include', '../external/boost_1_46_1'],
                                         #, ignore_gccxml_output = True
                                         #, indexing_suite_version = 1
                                         gccxml_path=gccxml_path,
                                         )

        print commands.getoutput('rm %s/exposed_decl.pypp.txt' % fname_base)

        depend_on = IncludeDict.get(fl, (None, None, []) )[2][:]
        #depend_on.extend( last_header_list )
        #if None in depend_on: depend_on = []
        print 'depend_on = %s' % depend_on
        if depend_on:
            depend_on = depend_on[-1:]
            decls_file = file(fname_base + '/exposed_decl.pypp.txt', 'w')
            for f in depend_on:
                sp_name = f.split('/')
                name = 'python/bindings/rosetta/' + '/'.join(sp_name[:-1]) + '/_' + sp_name[-1].replace('.', '_')[:-3]
                decls_file.write( file(name + '.exposed_decl.pypp.txt').read() )

            decls_file.close()
            mb.register_module_dependency(fname_base)


        #print 'Excluding stuff... ---------------------------------------------'
        exclude.exclude(path, mb, hfile=fl)
        #mb.build_code_creator( module_name = '_' + dname )
        #def extr(something): return '"ABCDEF"'

        #mb.build_code_creator( module_name = '_'+hbase, doc_extractor=extr)#, doc_extractor=doxygen.doxygen_doc_extractor() )
        mb.build_code_creator( module_name = '_'+hbase, doc_extractor=tools.DoxygenExtractorPyPP.doxygen_doc_extractor() )

        mb.code_creator.user_defined_directories.append( os.path.abspath('.') )  # make include relative
        mb.write_module( os.path.join( os.path.abspath('.'), fname ) )

        exclude.finalize(fname, dest, path, mb, module_name='_'+hbase, add_by_hand = (fl==headers[0]), files=[fl])
        finalize_init(fname)

        del mb  # to free extra memory before compiling

        # saving generated 'exposed_decl.pypp.txt' for possible later use...
        commands.getoutput('cp %s/exposed_decl.pypp.txt %s' % (fname_base, decl_name))

        #exclude.finalize_old(fname, path, mb)  # remove init for some reasons.
        #print 'Module name="%s"' % dname

        # Mac OS compiling options: -pipe -ffor-scope -W -Wall -pedantic -Wno-long-long
        #        -Wno-long-double -O3 -ffast-math -funroll-loops -finline-functions
        #        -finline-limit=20000 -s -Wno-unused-variable -march=prescott -fPIC

        if not Options.one_lib_file:
            if execute("Compiling...", # -fPIC
                "gcc %(fname)s -o %(obj_name)s -c \
                 %(add_option)s -I../external/include \
                 %(include_paths)s " % dict(add_option=add_option, fname=fname, obj_name=obj_name, include_paths=include_paths),
                 True):

                print 'Compiliation failed... Creating empty C++ file to test includes...'
                fname =    fname_base + '/' + '_0_.cc'
                obj_name = fname_base + '/' + '_0_.o'
                f = file(fname, 'w');  f.write('#include <%s>\n' % fl);  f.close()
                execute("Compiling empty C++...",
                    "gcc %(fname)s -o %(obj_name)s -c \
                     %(add_option)s -I../external/include \
                     %(include_paths)s " % dict(add_option=add_option, fname=fname, obj_name=obj_name, include_paths=include_paths))

                sys.exit(1)



#                 -lObjexxFCL -lutility -lstdc++ \
            #all_libs = '-lObjexxFCL -lutility -lnumeric -lcore -lprotocols'
            #if Platform == 'cygwin': all_libs = '-lmini'

#-lObjexxFCL -lutility -lnumeric -lcore -lprotocols -lstdc++ \
            #if Platform == 'cygwin': runtime_libpaths = ''
            execute("Linking...", # -fPIC -ffloat-store -ffor-scope
                "cd %(dest)s/../ && gcc %(obj)s %(add_option)s  \
                -lmini -lstdc++ \
                 -l%(python_lib)s \
                 -l%(boost_lib)s \
                 %(libpaths)s %(runtime_libpaths)s -o %(dst)s" % dict(add_option=add_loption, obj=obj_name,
                    dst=dst_name, libpaths=libpaths, runtime_libpaths=runtime_libpaths, dest=dest, boost_lib=Options.boost_lib,
                    python_lib=Options.python_lib,
                    )
                 )



    # Generate just one lib file -------------------------------------------------------------------
    if Options.one_lib_file:
        all_files = cc_files
        max_files = Options.one_lib_max_files or len(all_files)

        f = file( dest + '/' + path + '/__init__.py', 'w'); f.close()  # truncating the __init__ file

        for f_range in range(0, len(all_files), max_files):
            cc_files = all_files[ f_range : f_range+max_files ]
            #print 'cc_files:', cc_files


            dir_base = 'rosetta_%s_%03d' % (path.replace('/', '_'), f_range/max_files)
            cc_all   = fname_base + '/' + '_%s.cc'  % dir_base
            obj_name  = fname_base + '/' + '_%s.o'  % dir_base
            dst_name = fname_base + '/' + '_%s.so' % dir_base
            if Platform == 'cygwin' : dst_name = fname_base + '/' + '_%s.dll' % dir_base

            begining, end = '', ''
            for f in cc_files:
                lines = file(f).read().split('\n')
                b, e, = 9999999, 9999999

                for i in range( len(lines) ):
                    if lines[i].startswith('BOOST_PYTHON_MODULE('): b = i
                    if lines[i].startswith('}'): e = i

                begining += '\n'.join( lines[:b] )
                end += '\n'.join( lines[b+1:e] )

            text = begining + 'BOOST_PYTHON_MODULE(_'+ dir_base + ') {\n' + end + '\n}\n'

            f = file(cc_all, 'w');  f.write(text);  f.close()

            cc_recompile = True

            if Options.update:
                cc_recompile = False
                try:
                    for f in cc_files:
                        print 'checking file %s...' % f
                        if os.path.getmtime(f) >= os.path.getmtime(dst_name):
                            cc_recompile = True
                            break
                except os.error: cc_recompile = True

            if cc_recompile:
                execute("Compiling one lib...",
                    "gcc %(fname)s -o %(obj_name)s -c \
                     %(add_option)s -I../external/include \
                     %(include_paths)s " % dict(add_option=add_option, fname=cc_all, obj_name=obj_name, include_paths=include_paths) )

                #all_libs = '-lObjexxFCL -lutility -lnumeric -lcore -lprotocols'
                #if Platform == 'cygwin': all_libs = '-lmini'

                execute("Linking one lib...",
                    "cd %(dest)s/../ && gcc %(obj_name)s %(add_option)s  \
                     -lmini -lstdc++ \
                     -l%(python_lib)s \
                     -l%(boost_lib)s \
                     %(libpaths)s %(runtime_libpaths)s -o %(dst)s" % dict(add_option=add_loption, obj_name=obj_name,
                        dst=dst_name, libpaths=libpaths, runtime_libpaths=runtime_libpaths, dest=dest, boost_lib=Options.boost_lib,
                        python_lib=Options.python_lib,
                        )
                     )
            #print 'Finalizing Creating __init__.py file...'
            f = file( dest + '/' + path + '/__init__.py', 'a');
            f.write('from _%s import *\n' % dir_base);  f.close()


    print 'Done!'
    return new_headers


def buildModule_All(path, dest, include_paths, libpaths, runtime_libpaths, gccxml_path):
    ''' Non recursive build buinding for given dir name, and store them in dest.
        path - relative path to namespace
        dest - path to root file destination, actual dest will be dest + path
    '''
    print 'Processing', path
    dname = os.path.basename(path)

    # Creating list of headers
    headers = [os.path.join(path, d)
                for d in os.listdir(path)
                    if os.path.isfile( os.path.join(path, d) )
                        and d.endswith('.hh')
                        and not d.endswith('.fwd.hh')
                        ]

    for i in exclude.exclude_header_list:
        if i in headers:
            print "Excluding header:", i
            headers.remove(i)

    print headers

    fname_base = dest + '/' + path

    if not os.path.isdir(fname_base): os.makedirs(fname_base)

    dst_name = fname_base + '/' + '_' + dname + '.so'

    if not headers:  # if source files is empty then __init__.py should be empty too
        print 'Creating __init__.py file...'
        f = file( dest + '/' + path + '/__init__.py', 'w');
        f.close()
        return

    include_paths = ' -I'.join( [''] + include_paths + ['../src/platform/linux' , '../src'] )

    nsplit = 99999  # 4 files per iteration
    fcount = 0
    namespace_wraper = ' '.join( [ 'namespace '+ n + ' { ' for n in path.split('/')] ) + ' %s ' + ' '.join( [ '}' for n in path.split('/')] ) + '\n'

    while headers:
        files = headers[:nsplit]
        headers = headers[nsplit:]
        #print 'Binding:', files

        fname = fname_base +    '/' + '_' + dname + '.cc'
        inc_name = fname_base +    '/' + '_' + dname + '.hh'
        obj_name = fname_base + '/' + '_' + dname + '.o'

        # Read all the include files and write them in to a single one. (Walk around for Py++ bug)
        #incList = ['#include <'+x+'>\n' for x in files]

        lines = []
        for f in files:
            lines.extend( open(f).read().split('\n') )

        begining, middle, rest = '', '', ''

        for l in lines:
            if l.startswith('#ifndef'): l = '// Commented by BuildBindings.py ' + l
            if l.startswith('#endif'): l = '// Commented by BuildBindings.py ' + l
            if l.startswith('#define'):
                begining += l + '\n'
                l = '// Moved to top of the file ' + l

            if l.startswith('#include'):
                middle += l + '\n'
                l = '// Moved to middle of the file ' + l

            rest += l + '\n'

        S = re.findall(r'struct (.*){', rest)
        for s in S:
            begining += namespace_wraper % ('struct ' + s + ';\n' )

        S = re.findall(r'class (\w*)', rest)
        for s in S:
            begining += namespace_wraper % ('class ' + s + ';\n' )


        h = open(inc_name, 'w')
        h.write( begining + '\n' + middle + '\n'+ rest + '\n')
        h.close()

        mb = module_builder.module_builder_t(files=files #[inc_name]  #files=files
                                         , include_paths = ['../src/platform/linux', '../external/include', '../external/boost_1_46_1'],
                                         #, ignore_gccxml_output = True
                                         #, indexing_suite_version = 1
                                         gccxml_path=gccxml_path,
                                         )
        print 'Excluding stuff... ---------------------------------------------'
        exclude.exclude(path, mb)

        mb.build_code_creator( module_name = '_' + dname )

        mb.code_creator.user_defined_directories.append( os.path.abspath('.') )  # make include relative
        mb.write_module( os.path.join( os.path.abspath('.'), fname ) )

        exclude.finalize(fname, dest, path, mb)

        del mb  # to free extra memory before compiling

        #exclude.finalize_old(fname, path, mb)  # remove init for some reasons.
        print 'Module name="%s"' % dname

        if Platform == 'linux':
            add_option = '-malign-double' if PlatformBits == '32' else '-fPIC'

        execute("Compiling...", # -fPIC
            "gcc %(fname)s -o %(obj_name)s -c \
             %(add_option)s -ffloat-store -ffor-scope -I../external/include \
             %(include_paths)s " % dict(add_option=add_option, fname=fname, obj_name=obj_name, include_paths=include_paths) )

        fcount += 1

    obj_name = fname_base + '/' + '_' + dname + '.o '
    #for i in range(fcount):
    #    obj_name += fname_base + '/' + '_' + dname + '.' + str(i) + '.o '

    libpaths = ' -L'.join( ['', dest] + libpaths )
    runtime_libpaths = ' -Xlinker -rpath '.join( [''] + runtime_libpaths + ['rosetta'] )

    if Platform == 'linux': add_option = '-shared'
    else: add_option = '-dynamiclib'

    global Options
    execute("Linking...",
            "cd %(dest)s/../ && gcc %(add_option)s -fPIC \
             -ffloat-store -ffor-scope \
             -lObjexxFCL -lutility -lnumeric -lcore -lprotocols -lstdc++ \
             -l%(python_lib)s \
             -l%(boost_lib)s \
             %(libpaths)s %(runtime_libpaths)s %(obj)s -o %(dst)s" % dict(add_option=add_option, obj=obj_name,
                dst=dst_name, libpaths=libpaths, runtime_libpaths=runtime_libpaths, dest=dest, boost_lib=Options.boost_lib),
                python_lib=Options.python_lib,
             )

# boost_python-xgcc40-mt-1_37

#              %(libpaths)s %(runtime_libpaths)s %(obj)s -o %(dst)s.dylib" % dict(obj=obj_name,

#     -lboost_python-gcc40-mt-1_36


    print 'Done!'
#             -Xlinker -rpath %(runtime_libpaths)s -Xlinker -rpath /home/sergey/y/lib \
# -lboost_python-gcc34-mt-1_34_1




'''

print 'Checking for ParallelPython...',
try:
    import pp
    ParallelPython = True
    print ' found!'
    JobServer = pp.Server()
    Jobs = []
except ImportError:
    ParallelPython = False
    print ' not found!'


def build():
    mb = module_builder.module_builder_t(
         files=['src/utility/exit.hh'])

        #, gccxml_path=gccxml.executable ) #path to gccxml executable

    #mb.print_declarations()



    mb.build_code_creator( module_name='rosetta', doc_extractor=doxygen.doxygen_doc_extractor() )
    mb.code_creator.user_defined_directories.append( os.path.abspath('.') )  # make include relative

    mb.write_module( os.path.join( os.path.abspath('.'), 'generated.cpp' ) )


'''

if __name__ == "__main__": main(sys.argv)

# class revision 26929
# ? Score Function, Conformation?
