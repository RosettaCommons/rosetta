#!/bin/bash 
## Updated 22 Feb 2015 by Vikram K. Mulligan, Baker Lab, for latest
## Argonne Blue Gene/Q system "Mira" setup.
cd ../
./make_project.py all
cd build_staticmpi_bluegene_clang
#export CXX=mpicxx #For gcc compilation
export CXX=mpiclang++11 #For bgclang compilation
/soft/buildtools/cmake/current/bin/cmake ./
make -j 40

##For gcc compilation, make sure that your .soft file contains the following:
#    +mpiwrapper-gcc
#    +bgqtoolchain-gcc484
#    +python
# @default

##For bgclang compilation, the .soft file should include:
#      +mpiwrapper-bgclang
#      +python
# @default
