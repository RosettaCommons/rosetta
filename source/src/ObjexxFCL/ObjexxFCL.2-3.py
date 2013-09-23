#!/usr/bin/env python

# ObjexxFCL Version 2 to 3 Code Migration Tool
#
# Project: ObjexxFCL
#
# Version: 1.0.0.Rosetta
#
# Language: Python
#
# Author: Stuart G. Mentzer (Objexx Engineering, Inc.)
#
# Copyright (c) 2009 Objexx Engineering, Inc.


# Imports
import getopt, glob, os, re, sys


# Regex search and replace pairs
regexs = [
 ( re.compile( r'([^A-Za-z0-9_:]|^)CPArray([<_ ]|\.fwd\.hh|\.hh)' ), r'\1CArrayP\2' ),
 ( re.compile( r'([^A-Za-z0-9_:]|^)FArrayB([<_ ]|\.fwd\.hh|\.hh)' ), r'\1FArray\2' ),
 ( re.compile( r'([^A-Za-z0-9_:]|^)FArray([1-6])DB([<_ ]|\.fwd\.hh|\.hh)' ), r'\1FArray\2\3' ),
 ( re.compile( r'([^A-Za-z0-9_:]|^)FArray([1-6])Da([<_ ]|\.fwd\.hh|\.hh)' ), r'\1FArray\2A\3' ),
 ( re.compile( r'([^A-Za-z0-9_:]|^)FArray([1-6])Dp([<_ ]|\.fwd\.hh|\.hh)' ), r'\1FArray\2P\3' ),
 ( re.compile( r'([^A-Za-z0-9_:]|^)FArray([1-6])D(\.io\.hh)' ), r'\1FArray\2\3' ),
 ( re.compile( r'ObjexxFCL/internal' ), r'ObjexxFCL' ),
 ( re.compile( r'ObjexxFCL(/format|)/formatted\.(i|o|io)\.hh' ), r'ObjexxFCL/format.hh' ),
 ( re.compile( r'ObjexxFCL/FArray\.fwd\.hh' ), r'ObjexxFCL/FArray.all.fwd.hh' ),
 ( re.compile( r'ObjexxFCL/FArray([1-6])Ds(\.fwd\.hh|\.hh)' ), r'ObjexxFCL/FArray\1.all\2' ),
]


# Main
def main():

    # Get options and arguments
    try:
        opts, args = getopt.getopt( sys.argv[1:], 'ht', [ 'help', 'tree' ] )
    except getopt.GetoptError, err:
        print str( err )
        help_display()
        sys.exit( 2 )
    tree = False
    for o, a in opts:
        if o in ( '-h', '--help' ): # Help
            help_display()
            sys.exit()
        elif o in ( '-t', '--tree' ): # Walk the current source tree
            tree = True
            if len( args ) > 0:
                print 'ERROR: The --tree option cannot be combined with a C++ file argument'
                sys.exit( 2 )
        else: # Treat as the input file
            raise Exception, 'Unsupported option: ' + o

    # Process arguments
    if tree:
        cwd = os.getcwd()
        if ( cwd in ( '/', '\\' ) ) or ( ( os.name == 'nt' ) and ( len( cwd ) == 3 ) and ( cwd[ 1:3 ] == ':\\' ) ):
            print 'ERROR: Migrating the tree from a root directory is not allowed'
            sys.exit( 2 )
        for root, dirs, files in os.walk( cwd ):
            for file in files:
                base, ext = os.path.splitext( os.path.basename( file ) )
                if ext.startswith( '.' ): ext = ext[1:]
                if ext in ( 'cc', 'cpp', 'cxx', 'c++', 'C', 'h', 'hh', 'hpp', 'hxx', 'h++', 'i', 'ii', 'ipp', 'ixx', 'i++' ):
                    if root == '.':
                        iname = file
                    elif ( root[:2] == '.\\' ) or ( root[:2] == './' ):
                        iname = os.path.join( root[2:], file )
                    else:
                        iname = os.path.join( root, file )
                    process_file( iname )
    else:
        if len( args ) == 0:
            print 'ERROR: No C++ file(s) specified'
            print
            help_display()
            sys.exit( 2 )
        for arg in args:
            for iname in glob.glob( arg ):
                process_file( iname )


def process_file( iname ):

    # Open input file
    if not os.path.exists( iname ):
        for ext in ( 'cc', 'cpp', 'cxx', 'c++', 'C', 'h', 'hh', 'hpp', 'hxx', 'h++', 'i', 'ii', 'ipp', 'ixx', 'i++' ):
            if os.path.exists( iname + '.' + ext ):
                iname = iname + '.' + ext
                break
        print 'ERROR: No C++ file found matching the name ' + str( iname )
        sys.exit( 2 )
    base, ext = os.path.splitext( os.path.basename( iname ) )
    if ext.startswith( '.' ): ext = ext[1:]
    if ext not in ( 'cc', 'cpp', 'cxx', 'c++', 'C', 'h', 'hh', 'hpp', 'hxx', 'h++', 'i', 'ii', 'ipp', 'ixx', 'i++' ):
        print 'WARNING: C++ file does not have a standard extension: ' + str( iname )
    ifile = open( iname, 'r' )

    # Read and close input file
    ilines = ifile.readlines()
    ifile.close()

    # Migrate each line
    changed = False
    olines = []
    for iline in ilines:
        oline = iline
        for regex in regexs:
            oline = regex[0].sub( regex[1], oline )
        if oline != iline:
            changed = True
        olines.append( oline )

    # Write modified file if changed
    if changed: # Migrated
        print iname, ' migrated'
        ofile = open( iname, 'w' )
        ofile.writelines( olines )
        ofile.close()
    else: # No changes
        print iname, ' unaffected'


# Help Display
def help_display():
    print 'ObjexxFCL 2.x to 3.x Migration Tool'
    print
    print 'Renamed Types:'
    print ' CPArray   -> CArrayP'
    print ' FArrayB   -> FArray    (base class of all FArrays)'
    print ' FArrayNB  -> FArrayN   (base class of N-D FArrays)'
    print ' FArrayNDa -> FArrayNA  (argument FArrays)'
    print ' FArrayNDp -> FArrayNP  (proxy FArrays)'
    print '(FArrayND is still used for real FArrays)'
    print
    print 'Merged: formatted.i/o/io.hh/cc -> format.hh/cc'
    print
    print 'Usage:'
    print
    print os.path.basename( sys.argv[ 0 ] ), '[options] <C++_file(s)>'
    print
    print '  -h  --help   Help display'
    print
    print '  -t  --tree   Migrate all C++ files in the current tree'
    print
    print 'WARNING: THIS TOOL PERFORMS IN-PLACE FILE CONVERSION!!!'


# Runner
if __name__ == '__main__':
    main()
