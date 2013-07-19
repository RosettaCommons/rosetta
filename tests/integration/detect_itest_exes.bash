#!/bin/bash
#The purpose of this script is to delete integration tests from release packages.  It deletes tests that rely on unreleased code.  DO NOT RUN THIS SCRIPT EVER (except if you are preparing a release).  

for directory in `ls tests`
do

#this monster finds all of the lines containing executable names from the tests' command files (first grep), strips the line down to only the part where the name is plus the flanking bin and binext string replacement junk (egrep), then also removes the string replacement junk (sed commands), then sorts and uniqs the list
for executable in `grep "%(bin)s" tests/$directory/command | grep -v MY_MINI_PROGRAM | egrep -o '%\(bin\)s/.*%\(binext\)s' | sed 's/%(bin)s\///g' | sed 's/.%(binext)s//g' | sort | uniq`
do
#echo $directory $executable

    if grep -q $executable ../../source/src/apps.src.settings
	then
	echo $executable "found in apps.src.settings" $directory "may be ok"
	else
	echo $executable "not found, REMOVE" $directory
	git rm -r tests/$directory #--dry-run
	echo "rm -r ref/$directory"
    fi
done

done