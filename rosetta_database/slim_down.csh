#!/bin/bash 

for i in $( find . -name '*.slim' ); do 
  echo copy $i to $( echo $i | sed s/.slim// )
  cp $i $( echo $i | sed s/.slim// ) 
done
