#!/usr/bin/env bash
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## prun.sh: same as rrun.sh without the ppk starting structure to rosetta, so we can run from the native for prepacking

# usage:
# rosettarun.sh compiler pdb_code

# makes a directory named for the pdb_code and runs in that directory.
# works with a Condor script or alone.

if [ $# -lt 3 ]
then
  echo ERROR--Not enough arguments: $@
  echo Usage: $0 compiler pdb_code prefix_or_number extra_flags
  exit 1
fi

InputLine=$@
COMPILER=$1
PDB=$2
prefix=$3
shift 3

#--- error checking -------------
#if [ ! -f paths.txt ]; then echo No paths.txt ...exiting; exit;fi
if [ ! $prefix ]; then echo No prefix detected...exiting; exit;fi

#====================================================================

cd $PDB

#--- Executable -----------------
exe=$COMPILER
if [ ! -x $exe ]; then echo No executable $exe; exit; fi


#--- Arguments ------------------
# standard arguments plus optional arguments
args="$@"

echo --------------------------------
echo Rosetta Run in $PWD $(eval date)
echo $(eval uname -a)
echo $(basename $0) $InputLine
echo PDB: $PDB
echo prefix: $prefix
echo pdbpath: $pdbpath
echo exe: $exe
echo args: $args
echo --------------------------------


#--- Run ------------------------
fail="no"
if [ $CONDOR_VM ] && ! echo $args | grep -qE "\-verbose|\-inform|\-chat|\-yap"
then
    # Running inside condor, don't keep output unless -verbose flag is on
    /usr/bin/time nice $exe $args 2>/dev/null || fail=yes
else
    # Running interactively
    /usr/bin/time nice $exe $args 2>&1        || fail=yes
fi


echo -------------------
if [ "$fail" = "no" ]; then 
echo Run finished \($PDB\)
else
echo Run failed!! \($PDB\)
fi
echo -------------------


#--- Clean -----------------------
# tar and move files from scratch directories
#if echo $pdbpath | grep scratch
#then
#    savedir=$PWD
#    echo pdb files in $pdbpath are being zipped and moved to $savedir
#    cd $pdbpath
#    tarfile=$prefix$PDB.pdbs.tar.gz
#    tar -czvf $tarfile *$prefix*pdb --remove-files
#    mv $tarfile $savedir
#fi

cd ../

exit
