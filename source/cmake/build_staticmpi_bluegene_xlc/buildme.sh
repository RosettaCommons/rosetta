#!/bin/bash 
sed -i 's/python2\.7/python2/g' ../*.py ../../*.sh
cd ../
./make_project.py all
cd build_staticmpi_bluegene_xlc
export CXX=/bgsys/drivers//ppcfloor/comm/bin/xl/mpixlcxx
/soft/buildtools/cmake/current/bin/cmake ./
make -j 40
sed -i 's/python2/python2.7/g' ../*.py ../../*.sh
