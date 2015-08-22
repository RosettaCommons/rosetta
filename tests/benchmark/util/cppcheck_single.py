#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   cppcheck_single.py
## @brief  A helper script to update cppcheck cached files on one or more Rosetta files
##     Usage: "cppcheck_single.py "compile_type" filenames ...
##     Needs to be run from within the Rosetta/main/source/src directory
## @author Rocco Moretti (rmorettiase@gmail.com)

import os, os.path, commands
import sys

INCLUDES = "-I ./ -I ./platform/linux/ -isystem ../external/boost_1_55_0/ -isystem ../external/include/ -isystem ../external/dbio/"

CACHE_DIRECTORY = '../build/cppcheck/src/' # from the source/src/ directory

TESTS_TO_RUN = "warning,style,performance" # "warning,style,performance,portability" -- error is always on

IGNORED_DIRS = ['./devel','./apps/pilot', './python', './protocols/sparta']
IGNORED_FILES = []

def parse_defines( compile_type ):
    extras = compile_type.split('-')
    defines = ["BOOST_ERROR_CODE_HEADER_ONLY","BOOST_SYSTEM_NO_DEPRECATED","__USE_XOPEN2K8"]
    if 'cxx11' or 'cxx11thread' in extras:
        defines.append( "PTR_STD" )
    else:
        defines.append( "PTR_BOOST" )

    if 'cxx11' in extras:
        defines.extend( ["CXX11"] )
    if 'cxx11thread' in extras:
        defines.extend( ["CXX11", "MULTI_THREADED"] )
    if 'mpi' in extras:
        defines.extend( ["USEMPI"] )
    if 'mysql' in extras:
        defines.extend( ["USEMYSQL"] )
    if 'postgres' in extras:
        defines.extend( ["USEPOSTGRES"] )
    if 'graphics' in extras:
        defines.extend( ["GL_GRAPHICS"] )

    return '-D ' + ' -D '.join( defines )

def get_cache_filename(filename, compile_type ):
    return os.path.join( CACHE_DIRECTORY, compile_type, filename ) + '.cppcheck'

# Check if the cache file is up-to-date
def check_cache_file( filename, compile_type ):
    cache_filename = get_cache_filename(filename, compile_type)
    if not os.path.exists( cache_filename ):
        return

    last_run = os.path.getmtime( cache_filename )

    if os.path.getmtime( filename ) > last_run or os.path.getmtime( "cppcheck_suppressions.txt" ) > last_run:
        os.remove( cache_filename )
        return

    define_options = parse_defines( compile_type )
    commandline = 'clang++ -MM {0} {1} {2}'.format(filename, define_options, INCLUDES)

    res, output = commands.getstatusoutput( "{} 2> /dev/null".format( commandline ) )
    if res:
        os.remove( cache_filename )
        return

    output = output.split()
    for fn in output[1:]:
        if fn == '\\':
            continue
        if os.path.getmtime( fn ) > last_run:
            os.remove( cache_filename )
            return

#If the cached file is not up-to-date, then run cppcheck on the parent file.
def process_file( filename, compile_type ):
    cache_filename = get_cache_filename(filename,compile_type)
    if os.path.exists( cache_filename ):
        print "UP TO DATE:", filename
        return
    print "UPDATING:", filename
    tests_to_run = TESTS_TO_RUN

    define_options = parse_defines( compile_type )
    commandline = "cppcheck --suppressions cppcheck_suppressions.txt {defines} --enable={tests} {filename}".format(defines=define_options,tests=tests_to_run,filename=filename )

    res, output = commands.getstatusoutput( commandline )

    if res:
        selected_lines = "Error running cppcheck: " + commandline + "\n\n" + output
    else:
        selected_lines = [ line.strip() for line in output.splitlines() if line.startswith('[') ]

    #Filename may contain directories - we need to make sure they're present first.
    dirname = os.path.dirname( cache_filename )
    if not os.path.exists( dirname ):
        os.makedirs( dirname )
    with open( cache_filename, 'w' ) as f:
        if selected_lines:
            f.write( '\n'.join(selected_lines) + '\n' )

if __name__ == "__main__":
    compile_type = sys.argv[1]
    for filename in sys.argv[2:]:
        if filename in IGNORED_FILES:
            print "SKIPPING", filename
            continue
        skip = False
        for dirname in IGNORED_DIRS:
            if filename.startswith( dirname ):
                skip = True
                print "SKIPPING", filename
                break
        if not skip:
            print filename
            check_cache_file( filename, compile_type )
            process_file( filename, compile_type )
