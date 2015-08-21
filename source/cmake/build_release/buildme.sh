#!/bin/bash 
## This script will build Rosetta in release mode using gcc and 4 processes.
## If you want more parallel compilation processes, change the "make -j <number>" line.
cd ../
./make_project.py all
cd build_release/
cmake ./
make -j 4
