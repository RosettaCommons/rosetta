#!/usr/bin/env python

# gen_cmake_build.py - script to parse SCons configuration and generate
# requisite build files. You should be able to make this call to the script:
# gen_cmake_build.py gcc,macos,4.0 build_mac_gcc and get a CMake build
# directory in build_mac_gcc.

from __future__ import print_function

import re
import sys

if len( sys.argv) == 1:
	print("usage: gen_cmake_build.py build_tags")
	sys.exit(1)
tag = sys.argv[1]

build_tokens = sys.argv[1:len(sys.argv)]
build_tokens.sort()

# load basic.settings.
settings_file = '../tools/build/basic.settings'
str = ''.join( file( settings_file ).readlines() )
exec(str)

def remove_whitespace(s):
	return re.sub( r'\s', '', key )

for key in settings.keys():
	tokens = remove_whitespace(key).split(',')
	tokens.sort()
	if tokens == build_tokens:
		print(settings[key])
