#!/bin/bash 
cd ../
make_project.py all
cd build_staticmpi_bluegene
export CXX=/bgsys/drivers/ppcfloor/comm/bin/mpicxx
cmake ./
make -j 8
