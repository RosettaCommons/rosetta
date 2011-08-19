#!/bin/bash 

for i in $( find . -name '*.slim' ); do 
  svn revert $( echo $i | sed s/.slim// ) 
done
