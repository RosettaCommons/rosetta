#!/bin/bash

export CXX=`which clang++`
export CC=`which clang`
export CCC_CXX=`which clang++`
export CCC_CC=`which clang`

pushd ../
./make_project.py all
popd
scan-build cmake -G "Unix Makefiles" .
mkdir -p ../../build/clang_SA/
scan-build --keep-going -plist-html --keep-empty -o ../../build/clang_SA/ make -j ${1}

# We DON'T use the --status-bugs option for scan-build, so the exit code will be
# the exit code of the build system
exit

