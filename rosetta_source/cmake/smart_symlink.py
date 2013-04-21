#!/usr/bin/python

# smart_symlink.py - quick script for symlinking binaries into mini/bin
# usage: python ../smart_symlink.py gcc release

import sys
import glob
import os

platform = sys.platform

if platform == 'darwin':
    platform = 'macos'
elif platform.startswith('linux'):
    platform = 'linux'

compiler = sys.argv[1]
mode     = sys.argv[2]
files    = sys.argv[3:]
path     = '../../bin/'
binext   = platform + compiler + mode

#files = glob.glob( '../build/apps/*.cmake' )

if not os.path.exists( path ):
  os.mkdir( path )

for filename in files:
    cmake_name = os.path.basename(filename)
    executable_name = cmake_name.split('.')[0]
    new_name = path + executable_name + '.' + binext
    default_name = path + executable_name + '.default.' + binext
    default_name_no_extension = path + executable_name

    if os.path.exists( executable_name ):
        # removing this print statement as it is unnecessary, and adds extra lines to compact ninja output.
        # print 'about to try symlinking ', executable_name, ' to ', new_name
        for new_name_for_symlink in [ new_name, default_name, default_name_no_extension ]:
            if os.path.exists( new_name_for_symlink ):
                os.remove( new_name_for_symlink )
                # print "unlinking",new_name_for_symlink
            #print 'symlinking ', executable_name, ' to ', new_name_for_symlink
            #os.symlink( os.path.abspath( executable_name ), new_name_for_symlink )
    else:
         print "Warning: %s doesn't exist!" % ( executable_name )
