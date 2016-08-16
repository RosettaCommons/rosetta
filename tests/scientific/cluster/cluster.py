#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   features.py
## @brief  A script to help debug/run cluster scientic tests
## @author Sergey Lyskov


import sys, commands, subprocess, re

import shutil
from os import path, environ, listdir, makedirs, symlink, readlink
from optparse import OptionParser, IndentedHelpFormatter


def main(argv):
    '''
A script to help debug/run cluster scientic tests

Each test has its own subdirectory in tests/, which serves as its name.
Each test has a files named "submit" and "analyze", which should have the command to be executed.
Variable substitution is possible using Python printf format; see examples.
To create a new test, just make a new directory with a "command" file and other needed files.
See tests/HOW_TO_MAKE_TESTS for more information.

Intended usage is to run the tests once to submit jobs to cluster: './cluster submit docking'
and when condor run is finished to run again to analyze results: './cluster analyze docking'

Please note that this script is for debuging/testing purposes only, it is not used in ScientificCluster daemon.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)
    parser.add_option("-d", "--database",
      default="", # processed below
      help="Path to Rosetta database. (default: $ROSETTA3_DB, ../../../database)",
    )

    parser.add_option("-m", "--mini_home",
      #default=path.join( path.expanduser("~"), "mini"),
      default = path.dirname( path.dirname( path.dirname( path.dirname(path.abspath(sys.argv[0])) ) ) ) + "/source",
      help="Directory where Rosetta is found (default: ../../../source)",
    )

    parser.add_option("--mode",
      default="release",
      help="In selecting binaries, which mode was used? (default: release)",
    )

    parser.add_option("-c", "--compiler",
      default="gcc",
      help="In selecting binaries, which compiler was used? (default: gcc)",
    )

    parser.add_option("--extras",
      default="default",
      dest="extras",
      help="in selecting binaries, which options were specified? (default: default)",
    )

    parser.add_option("--score-function",
      default="score12prime",
      dest="scfxn",
      help="Apply a score function parametrization in score_functions/ to the benchmark (default: score12prime)",
    )
    parser.add_option("--daemon", action="store_true", dest="daemon", default=False,
      help="generate daemon friendly output (off by default)"
    )

    parser.add_option("--lsf_queue_name",
      default="128cpu", action="store",
      help="The name of the queue when submitting jobs on the lsf cluster, to see available queues run 'bqueues' (use with 'submit_lsf' action). [default %default]"
    )

    parser.add_option("--num_cores",
      default=128, action="store",
      help="How many cores to request when submitting jobs. [default %default]"
    )

    parser.add_option("--output_dir",
      default="output", action="store",
      help="Base directory of where to deposit data and results. [default %default]"
    )

    parser.add_option("--run_type",
      default="condor", action="store",
      help="Indicate which type of the action should be run. Eg, 'lsf' means execute benchmark on a Load Sharing Facility cluster, while 'dryrun' means simulate the output of running the benchmark for debugging purposes. [default %default]"
    )

    (options, args) = parser.parse_args(args=argv)

    if options.database == parser.get_default_values().database:
        if environ.get('ROSETTA3_DB') is not None and \
                path.isdir(environ.get('ROSETTA3_DB')):
            options.database = environ.get('ROSETTA3_DB')
        elif path.isdir( path.join( path.expanduser("~"), "rosetta_database") ):
            options.database = path.join( path.expanduser("~"), "rosetta_database")
    	elif path.isdir("../../../database"):
            options.database = "../../../database"
        else:
            print "Can't find database at %s; please set $ROSETTA3_DB or use -d" % options.database
            return 1
    # Normalize path before we change directories!
    options.database = path.abspath(options.database)


    if len(args) != 2:
        print 'You must supplies action and test name in command line! For example: "./cluster.py submit docking" or ""./cluster.py analyze docking""'
        return

    valid_actions = ["submit", "analyze"]

    action, test = args

    if action not in valid_actions:
        print "ERROR: Action must be one of ['%s']" % "', '".join(valid_actions)
        return 1

    print 'Perform %s on test %s...' % (action, test)


    workdir = path.abspath( test )
    minidir = options.mini_home

    if not path.isdir(path.join("../score_functions", options.scfxn)) or \
            not path.isfile(path.join("../score_functions", options.scfxn, "flags")) or \
            not path.isfile(path.join("../score_functions", options.scfxn, "weights.wts")):
        print "ERROR: The specified score function '%s' must correspond to the directory:" % (options.scfxn)
        print "ERROR:   'tests/scientific/score_functions/%s'" % (options.scfxn)
        print "ERROR: containing the following files:"
        print "ERROR:   'flags' : a flags file specifying command line parameters"
        print "ERROR:   'weights.wts' : a score function weights file"
        print "ERROR: Available score functions are:"
        for scfxn_dir in sorted(listdir("../score_functions")):
            print "ERROR:   %s" % scfxn_dir
        return 1
    else:
        if path.isdir( path.join(workdir, test, options.scfxn) ):
            shutil.rmtree( path.join(workdir, test, options.scfxn) )
        copytree( path.join("../score_functions", options.scfxn), path.join(workdir, test, options.scfxn),
                  accept=lambda src, dst: path.basename(src) != '.svn' )


    database = options.database
    print 'Database:', database
    bin = path.join(minidir, "bin")
    pyapps = path.join(minidir, "src", "python", "apps")
    if sys.platform.startswith("linux"): platform = "linux" # can be linux1, linux2, etc
    elif sys.platform == "darwin": platform = "macos"
    else: platform = "_unknown_"
    extras = options.extras
    compiler = options.compiler
    mode = options.mode
    binext = extras+"."+platform+compiler+mode
    lsf_queue_name = options.lsf_queue_name
    num_cores = options.num_cores
    output_dir = options.output_dir
    run_type = options.run_type
    scfxn = options.scfxn
    # Read the command from the file "command"

    cmd = file(path.join(workdir, action)).read().strip()
    # cmd = cmd % vars() # variable substitution using Python printf style
    mvars = dict(minidir=minidir, database=database, workdir=workdir, platform=platform, bin=bin, compiler=compiler, mode=mode, binext=binext, lsf_queue_name=lsf_queue_name, num_cores=num_cores, output_dir=output_dir, run_type=run_type, scfxn=scfxn)
    cmd = cmd % mvars

    # Writing result to .sh file for future reference.
    f = file(path.join(workdir, action)+'.sh', 'w');  f.write(cmd);  f.close()

    # Creating python style file with arguments in case test want to use Python as script language
    f = file(path.join(workdir, '_arguments.py'), 'w');  f.write( str(mvars) );  f.close()

    # Executing created bash script...
    print commands.getoutput('cd %s && sh %s.sh' % (test, action) )


def copytree(src, dst, symlinks=False, accept=lambda srcname, dstname: True):
    """Recursively copy a directory tree using copy2(), with filtering.
    Copied from shutil so I could filter out .svn entries.
    """
    names = listdir(src)
    makedirs(dst)
    errors = []
    for name in names:
        srcname = path.join(src, name)
        dstname = path.join(dst, name)
        if not accept(srcname, dstname): continue
        try:
            if symlinks and path.islink(srcname):
                linkto = readlink(srcname)
                symlink(linkto, dstname)
            elif path.isdir(srcname):
                copytree(srcname, dstname, symlinks, accept)
            else:
                shutil.copy2(srcname, dstname)
            # XXX What about devices, sockets etc.?
        except (IOError, os.error), why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except shutil.Error, err:
            errors.extend(err.args[0])
    try:
        shutil.copystat(src, dst)
    except WindowsError:
        # can't copy file access times on Windows
        pass
    except OSError, why:
        errors.extend((src, dst, str(why)))
    if errors:
        raise shutil.Error, errors



if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


