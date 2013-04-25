#!/bin/bash 
cd ../
./make_project.py all
cd build_elscripts
export CC=/home/kenjung/bin/mpicc; export CXX=/home/kenjung/bin/mpicxx
cmake .
make -j 20
