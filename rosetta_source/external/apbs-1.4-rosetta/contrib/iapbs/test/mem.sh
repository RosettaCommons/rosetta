#!/bin/sh
#
# a quick check for memory leak
#
# usage:
# (../src/wrapper mol1-manual-loop.in &) ; ./mem.sh 

PID=`/usr/bin/pgrep -x wrapper`

LOG="memusage.txt.${PID}"

echo "# ElapsedTime VmSize VmRSS" > $LOG
echo "# " `date` >> $LOG

# seconds
PERIOD=0.1
counter=0.0

status_file=/proc/${PID}/status

echo "wrapper PID: $PID"

while [ -e $status_file ] ; do
   VM_SIZE=`awk '/VmSize/ {print $2}' < $status_file`
   VM_RSS=`awk '/VmRSS/ {print $2}' < $status_file`
   echo "$counter $VM_SIZE $VM_RSS" >> $LOG
   sleep $PERIOD
   counter=`echo $counter + $PERIOD|bc`
   VM_SIZE=""
   VM_RSS=""
done

echo "$0 script finished ..."

/usr/bin/gnuplot -persist <<EOF
#set xdata time
#set timefmt "%H:%M:%S"
set xlabel "time (s)"
set ylabel "memory (kB)"
plot "$LOG" using 1:2 w l t "vm size" , \\
 "$LOG" using 1:3 w l t "vm rss"
EOF

