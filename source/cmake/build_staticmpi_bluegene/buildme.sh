#!/bin/bash 
## Updated 22 Feb 2015 by Vikram K. Mulligan, Baker Lab, for latest
## Argonne Blue Gene/Q system "Mira" setup.
cd ../
./make_project.py all
cd build_staticmpi_bluegene
export CXX=/bgsys/drivers/ppcfloor/comm/gcc/bin/mpicxx
/soft/buildtools/cmake/current/bin/cmake ./
make -j 20

