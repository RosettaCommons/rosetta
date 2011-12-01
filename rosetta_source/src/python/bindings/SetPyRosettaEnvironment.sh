#!/bin/bash

# This script intended to set Pyrosetta environment variables so user can execute 'import rosetta' from
# any file system location. Use 'source SetPyRosettaEnvironment.sh' before starting to work with PyRosetta.

OLD_PATH=`pwd`

cd ${BASH_SOURCE[0]%/*}
PYROSETTA=`pwd`

#echo "Setting PyRosetta root as:" $PYROSETTA

export PYROSETTA
export PYTHONPATH=$PYROSETTA:$PYTHONPATH
export DYLD_LIBRARY_PATH=$PYROSETTA:$PYROSETTA/rosetta:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=$PYROSETTA/rosetta:$LD_LIBRARY_PATH
export PYROSETTA_DATABASE=$PYROSETTA/rosetta_database

cd $OLD_PATH
