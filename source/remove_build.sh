#!/bin/sh

# Run this script by providing the path to the build you want to remove,
# For instance to remove the old default build files I would type:
# "remove_build.sh build/src/release/linux/2.6/64/x86/gcc/4.5"

rm $1/*.so
rm $1/*default*
rm -rf $1/apps
rm -rf $1/basic
rm -rf $1/core
rm -rf $1/devel
rm -rf $1/numeric
rm -rf $1/ObjexxFCL
rm -rf $1/protocols
rm -rf $1/utility

