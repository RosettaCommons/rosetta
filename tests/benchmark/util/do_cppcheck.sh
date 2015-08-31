#!/bin/bash

# Assumes you're running it from within the Rosetta/main/source/src directory
# Need to set JOBS, and COMPILETYPE

JOBS=1
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

#Delete any cppcheck cached files from source files which no longer exist.
for f in `find ${CACHEDIR} -name '*.cppcheck'`; do
    ccfile=${f#${CACHEDIR}}
    ccfile=${ccfile%.cppcheck}
    if [ ! -e $ccfile ]; then
        rm $f
    fi
done

find ./ -name '*.cc' | sed "s|^|${CPPCHECK_DIR}/cppcheck_single.py ${COMPILETYPE} |g" > ${CACHEDIR}/commands.txt

../../tests/benchmark/util/parallel.py -q -j ${JOBS} ${CACHEDIR}/commands.txt > /dev/null # only the error output.

find ${CACHEDIR}/ -name '*.cppcheck' -exec cat {} \; > ${CACHEDIR}/all_lines.txt

#Running error lines, if any
grep -L '^\[' ${CACHEDIR}/all_lines.txt | uniq > ${CACHEDIR}/error_output.txt
#The problem list
grep '^\[' ${CACHEDIR}/all_lines.txt | sort | uniq > ${CACHEDIR}/output.txt

../../tests/benchmark/util/extract_lines.py ${CACHEDIR}/output.txt ../../tests/benchmark/util/cppcheck_known_lines.txt ${CACHEDIR}/oldissues.txt ${CACHEDIR}/newissues.txt

if [ -s ${CACHEDIR}/error_output.txt ]; then
    #error file size is zero, i.e. we have errors
    echo "ERRORS RUNNING CPPCHECK:"
    echo
    cat ${CACHEDIR}/error_output.txt
    echo
    echo "New Issues found:"
    echo
    cat ${CACHEDIR}/newissues.txt
    echo
    echo "Remaining historical issues:"
    cat ${CACHEDIR}/oldissues.txt
    echo
    echo "New issues found:" `cat ${CACHEDIR}/newissues.txt | wc -l`
    echo "Total issues found:" `cat ${CACHEDIR}/output.txt | wc -l`
    echo
    echo "Run walltime:" $(( $(date '+%s') - $starttime )) "s."
    exit 127
elif [ -s ${CACHEDIR}/newissues.txt ]; then
    #file size is non-zero, i.e. we have issues
    echo "New Issues found:"
    echo
    cat ${CACHEDIR}/newissues.txt
    echo
    echo "Remaining historical issues:"
    cat ${CACHEDIR}/oldissues.txt
    echo
    echo "New issues found:" `cat ${CACHEDIR}/newissues.txt | wc -l`
    echo "Total issues found:" `cat ${CACHEDIR}/output.txt | wc -l`
    echo
    echo "Run walltime:" $(( $(date '+%s') - $starttime )) "s."
    exit 1
elif [ -s ${CACHEDIR}/oldissues.txt ]; then
    #We don't have any new issues, but we do have old issues
    echo "Remaining historical issues:"
    cat ${CACHEDIR}/oldissues.txt
    echo
    echo "New issues found: None"
    echo "Total issues found:" `cat ${CACHEDIR}/output.txt | wc -l`
    echo
    echo "Run walltime:" $(( $(date '+%s') - $starttime )) "s."
    exit 0
else
    # No errors, success!
    echo "New issues found: None"
    echo "Total issues found: None"
    echo
    echo "Run walltime:" $(( $(date '+%s') - $starttime )) "s."
    exit 0
fi
