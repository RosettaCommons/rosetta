#!/bin/bash

echo 'Updating Rosetta options system...'

./update_options.sh

echo 'Updating ResidueType enum files...'

./update_ResidueType_enum_files.sh

echo 'Building documentation for each src/doxygen.*...'

mkdir html
cp src/index.html html/

rm -rf html/core+protocols html/all_else

for file in src/doxyfile_devel.* ; do
    echo "Building doc: $file"
    doxygen $file
done
