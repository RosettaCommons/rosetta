#!/bin/bash

# Get the directory that this script is in.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -z "$1" ]; then
   PREFIX="mpi_log"
else
   PREFIX=$1
fi

for log in ${PREFIX}*; do
    #echo "Cleaning " ${log}
    egrep -vf ${DIR}/ignore_list ${log} > ${log}.filt
    rm ${log}
done
