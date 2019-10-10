#!/bin/bash

set -e

export CXX=`which clang++`
export CC=`which clang`
export CCC_CXX=`which clang++`
export CCC_CC=`which clang`

pushd ../
./make_project.py all
popd
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON 
