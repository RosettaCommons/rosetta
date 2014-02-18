#!/usr/bin/env sh

doxygen_path=$(pwd)/$(dirname $0)
cd $doxygen_path

# Clean out the existing docs.
rm -rf html

# Build the new docs.
doxygen