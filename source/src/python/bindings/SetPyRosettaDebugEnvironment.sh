#!/bin/bash

# This script intended to set Pyrosetta environment variables so user can execute 'import rosetta' from
# any file system location. Use 'source SetPyRosettaEnvironment.sh' before starting to work with PyRosetta.

if [[ "${BASH_SOURCE[0]}" == "" ]]; then
    #echo "zsh like shell..."
    OLD_PATH=`pwd`
    PYROSETTA="$( cd "$( dirname "$0" )" && pwd )"/debug
else
    #echo "bash like shell..."
    OLD_PATH=`pwd`
    PYROSETTA="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/debug
fi


#if [[ "$0" == "bash"  ||  "$0" == "-bash" ]]; then
#    #echo "bash like shell..."
#    OLD_PATH=`pwd`
#    PYROSETTA="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#else
#    #echo "zsh like shell..."
#    OLD_PATH=`pwd`
#    PYROSETTA="$( cd "$( dirname "$0" )" && pwd )"
#fi

cd $PYROSETTA

#echo "Setting PyRosetta root as:" $PYROSETTA

export PYROSETTA
export PYTHONPATH=$PYROSETTA${PYTHONPATH+:$PYTHONPATH}
export DYLD_LIBRARY_PATH=$PYROSETTA:$PYROSETTA/rosetta${DYLD_LIBRARY_PATH+:$DYLD_LIBRARY_PATH}
export LD_LIBRARY_PATH=$PYROSETTA/rosetta${LD_LIBRARY_PATH+:$LD_LIBRARY_PATH}
export PYROSETTA_DATABASE=$PYROSETTA/rosetta_database

cd $OLD_PATH
