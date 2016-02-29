#!/bin/bash

RELEASE_MODE='devel'

[ "$1" == "--public" ] && RELEASE_MODE='release'

echo 'Updating Rosetta options system...'

./update_options.sh

echo 'Updating ResidueType enum files...'

./update_ResidueType_enum_files.sh

echo "Building documentation for each src/doxyfile_${RELEASE_MODE}.*..."

mkdir html
cp src/index.html html/

rm -rf html/core+protocols html/all_else

for file in src/doxyfile_${RELEASE_MODE}.* ; do
    echo "Building doc: $file"
    doxygen $file
done

