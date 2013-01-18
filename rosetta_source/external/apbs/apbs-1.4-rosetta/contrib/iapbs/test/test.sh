#!/bin/sh
#
# Run iAPBS validation tests.
#

wrapper=../src/wrapper
time="/usr/bin/time -v"
timing=./timing.log
touch $timing

export MCSH_HOME=/dev/null

# error for result comparison
ABSERR=1.0e-6

if [ "$1" = "single" ] ; then
    $wrapper apbs.in > apbs.out
    diff save/apbs.out.save apbs.out
    exit
fi

files="apbs apbs.d9 mol1-auto mol1-manual mol1-manual-loop \
 smpbe-ion smpbe-2ala apbs-forces apbs-forces-tot"

for i in $files
do
  echo "Working on $i ..."
  echo "Timing $i ..." 2>> $timing 1>&2
  $time $wrapper ${i}.in > ${i}.out 2>> $timing
  echo -n "Diffing ... "
  tmpfile=`mktemp ./tmp.XXX` || exit 1
#  diff save/${i}.out.save ${i}.out > $tmpfile
  awk -f ./ndiff.awk -v ABSERR=${ABSERR} \
      save/${i}.out.save ${i}.out > $tmpfile
  if [ -s $tmpfile ] ; then
      mv $tmpfile ${i}.out.diff
      echo "FAILED, see ${i}.out.diff."
  else
      echo "passed."
      rm $tmpfile
  fi
done

