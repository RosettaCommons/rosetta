#!/bin/bash 
cd ../
make_project.py all
cd build_staticmpi
cmake ./
make -j 20
