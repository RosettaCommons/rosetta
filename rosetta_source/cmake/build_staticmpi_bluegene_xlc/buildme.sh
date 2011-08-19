#!/bin/bash 
cd ../
make_project.py all
cd build_staticmpi_bluegene
export CXX=/bgsys/drivers/ppcfloor/comm/bin/mpixlcxx
cmake ./
make -j 8
