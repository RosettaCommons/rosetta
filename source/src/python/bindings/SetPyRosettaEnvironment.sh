#!/bin/bash

# This script intended to set Pyrosetta environment variables so user can
# execute 'import rosetta' from any file system location, even allowing
# the script to be a symbolic link.
#
# Use 'source SetPyRosettaEnvironment.sh' before starting to work with PyRosetta.

READLINK=""
if [ "Darwin" != $(uname -s) ]; then
    # FIXME: MacOS X does not have the -f option for readlink and realpath 
    #        is not available, either. 
    READLINK=$(which readlink)
fi

if [[ "${BASH_SOURCE[0]}" == "" ]]; then
    #echo "zsh like shell..."
    OLD_PATH=`pwd`
    if [ -n "$READLINK" ]; then
        PYROSETTA="$( dirname $( $READLINK -f "$0" ) )"
    else
        PYROSETTA="$( cd "$( dirname "$0" )" && pwd )"
    fi
else
    #echo "bash like shell..."
    OLD_PATH=`pwd`
    if [ -n "$READLINK" ]; then
        PYROSETTA="$( dirname $( $READLINK -f "${BASH_SOURCE[0]}" ) )"
    else
        PYROSETTA="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    fi
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

echo "Setting PyRosetta root as:" $PYROSETTA

export PYROSETTA
export PYTHONPATH=$PYROSETTA${PYTHONPATH+:$PYTHONPATH}
export DYLD_LIBRARY_PATH=$PYROSETTA:$PYROSETTA/rosetta${DYLD_LIBRARY_PATH+:$DYLD_LIBRARY_PATH}
export LD_LIBRARY_PATH=$PYROSETTA:$PYROSETTA/rosetta${LD_LIBRARY_PATH+:$LD_LIBRARY_PATH}
export PYROSETTA_DATABASE=$PYROSETTA/database

echo "Aliasing PyRosetta Toolkit GUI to pyrosetta_toolkit"
alias pyrosetta_toolkit='python $PYROSETTA/app/pyrosetta_toolkit/pyrosetta_toolkit.py'

cd -
