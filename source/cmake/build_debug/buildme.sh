#!/bin/bash 
## This script will build Rosetta in debug mode using gcc and 4 processes.
## If you want more parallel compilation processes, change the "make -j <number>" line.
cd ../
./make_project.py all
cd build_debug/
cmake ./
make -j 4

