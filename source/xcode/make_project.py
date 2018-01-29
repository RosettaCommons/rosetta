#!/usr/bin/env python



PATH_TO_ROOT = '../../'

import os, sys

sys.path.append( '..' )

# generate new version file
os.system( 'cd ..; python version.py' )

starting_directory = os.path.basename( os.getcwd() )
os.chdir( '..' )

os.chdir( starting_directory )

# (re)generate options files
os.system( 'cd ..; ./update_options.sh' )

# (re)generate ResidueType enum files
os.system( 'cd ..; ./update_ResidueType_enum_files.sh' )

# it's cmake time
os.system( 'cd ../cmake ; ./make_project.py all ; cd build_xcode ; cmake -G Xcode . ')

# it's symlink time
os.system( 'rm Rosetta.xcodeproj ; ln -s ../cmake/build_xcode/minirosetta.xcodeproj Rosetta.xcodeproj' ); 
