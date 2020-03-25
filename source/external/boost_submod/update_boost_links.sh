#!/bin/bash

# This script is necessary to update the boost/ subdirectory 
# What it will do is go through all of the submodules in this directory,
# and will symlink things into the directory.
# (It assumes that all of the submodules are updated appropriately.)

# Go to the boost directory
cd `dirname "${BASH_SOURCE[0]}"`

# Make sure we have the boost subdir
rm -rf boost
mkdir -p boost
cd boost

# Now link all the files in the submodule directories.
for f in ../*/include/boost/*; do
	bn=`basename $f`
	# These unfortunately are shared by multiple submodules, so have to be handled separately.
	if [ $bn == "detail" ]; then continue; fi
	if [ $bn == "utility" ]; then continue; fi
	if [ $bn == "exception" ]; then continue; fi
	if [ $bn == "functional" ]; then continue; fi
	if [ $bn == "pending" ]; then continue; fi
	if [ $bn == "numeric" ]; then continue; fi
	#echo "UPDATING: " ln -s -T ${f} ${bn}
	ln -s -T ${f} ${bn} # -T means don't interpret the second argument as a directory
done

# Combine the subdirs shared by several modules
for dir in detail utility exception functional pending numeric; do
    mkdir -p $dir
    for f in ../*/include/boost/${dir}/*; do
	bn=`basename $f`
	if [ ${dir}/${bn} == "pending/detail" ]; then continue; fi
	#echo "UPDATING: " ln -s -T ../${f} ${dir}/${bn}
	ln -s -T  ../${f} ${dir}/${bn}
    done
done

# Third level
for dir in "pending/detail"; do
   mkdir -p $dir
    for f in ../*/include/boost/${dir}/*; do
	bn=`basename $f`
	#echo "UPDATING: " ln -s -T ../../${f} ${dir}/${bn}
	ln -s -T  ../../${f} ${dir}/${bn}
    done
done

