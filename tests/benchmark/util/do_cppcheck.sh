#!/bin/bash

# Assumes you're running it from within the Rosetta/main/source/src directory
# Need to set JOBS, and COMPILETYPE

JOBS=16
COMPILETYPE=default

while getopts "j:e:" opt; do
  case $opt in 
    j)
      JOBS=$OPTARG
      ;;
    e)
      #COMPILETYPE assumes dash separated
      if ! [ -z "$OPTARG" ]; then
     	COMPILETYPE=$(echo $OPTARG | tr ',' '-')
      fi
      ;;
    '?')
      echo "Invalid option:-$OPTARG" >&2
      exit 127
      ;;
  esac
done

starttime=$(date '+%s')

# Check if cppcheck exists 
if ! which cppcheck &> /dev/null; then
    echo 'Cannot find cppcheck!'
    exit 127
fi

cppcheck_version=$(cppcheck --version)
if [ $? -ne 0 ]; then
    echo 'Problem running cppcheck! - ' ${cppcheck_version} 
    exit 127
fi 

echo "Running cppcheck tests with cppcheck version: " $cppcheck_version

CPPCHECK_DIR=../../tests/benchmark/util/
CACHEDIR=../build/cppcheck/src/${COMPILETYPE}/
mkdir -p ${CACHEDIR}

#find ./ -name '*.cc' 
find ./ -name '*.cc' | sed s'/^/\$\{CPPCHECK_DIR\}\/cppcheck_single.py \$\{COMPILETYPE\}/'g > ${CACHEDIR}/commands.txt

../../tests/benchmark/util/parallel.py -q -j ${JOBS} ${CACHEDIR}/commands.txt 

find ${CACHEDIR}/ -name '*.cppcheck' -exec cat {} \; > ${CACHEDIR}/all_lines.txt

#Running error lines, if any
grep -L '^\[' ${CACHEDIR}/all_lines.txt | uniq > ${CACHEDIR}/error_output.txt
#The problem list
grep '^\[' ${CACHEDIR}/all_lines.txt | sort | uniq > ${CACHEDIR}/output.txt

echo "Run walltime:" $(( $(date '+%s') - $starttime )) "s."
echo

if ! [ -s ${CACHEDIR}/error_output.txt ]; then
    #error file size is zero, i.e. we have errors
    cat ${CACHEDIR}/error_output.txt
    echo
    cat ${CACHEDIR}/output.txt
    echo
    echo "Issues found:" `cat ${CACHEDIR}/output.txt | wc -l`
    exit 127
elif ! [ -s ${CACHEDIR}/output.txt ]; then
    #file size is non-zero, i.e. we have issues
    cat ${CACHEDIR}/output.txt
    echo
    echo "Issues found:" `cat ${CACHEDIR}/output.txt | wc -l`
    exit 1
else
    # No errors, success!
    echo "Issues found: None"
    exit 0
fi
