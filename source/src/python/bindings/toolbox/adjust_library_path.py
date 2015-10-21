#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
# This code made available under dual license: RosettaCommons license and GPLv3

## @file   adjust_library_path.py
## @brief  Adjust shared library pasth in namespace .so files and set path it absolute value (required in El Capitan)
## @author Sergey Lyskov

import os, sys, os.path, subprocess

if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
elif sys.platform == "win32" : Platform = "windows"
else: Platform = "_unknown_"


def execute(message, command_line, return_=False, untilSuccesses=False, print_output=True, verbose=True):
    if verbose:
        print message
        print command_line

    while True:
        #print 'Platform:', Platform
        if Platform == 'cygwin':
            (res, output) = commands.getstatusoutput(command_line)
            if print_output or res: print output
            #sys.stdout.flush()
            #sys.stderr.flush()

        elif Platform == 'windows':
            # res = 0
            # try: output = subprocess.check_output(command_line, shell=True)
            # except subprocess.CalledProcessError as e: res, output = e.returncode, e.output
            # if print_output: print output

            po = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=sys.stdout, stderr=sys.stderr)
            while po.returncode is None: po.wait()
            res = po.returncode

        else:
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


        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if return_ == 'tuple': return(res, output)

    if res:
        if print_output: print "\nEncounter error while executing: " + command_line
        if not return_: sys.exit(1)

    if return_ == 'output': return output
    else: return res



for dir_name, _, files in os.walk('./rosetta'):
    #libmini = os.path.abspath('./rosetta/libmini.so')
    for f in files:
        if f.endswith('.so'):
            #print f, os.remove(dir_name+'/'+f)
            execute('', 'install_name_tool -change libboost_python.dylib `pwd`/libboost_python.dylib {}/{}'.format(dir_name, f), verbose=False)
            execute('', 'install_name_tool -change rosetta/libmini.dylib `pwd`/rosetta/libmini.dylib {}/{}'.format(dir_name, f), verbose=False)
