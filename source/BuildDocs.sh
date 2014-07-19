#!/bin/bash

echo 'Updating Rosetta options system...'

./update_options.sh

echo 'Populating ResidueType properties from database...'

./update_residue_properties.sh

echo 'Building documentation for each src/doxygen.*...'

mkdir html
cp src/index.html html/

rm -rf html/core+protocols html/all_else

for file in src/Doxyfile.* ; do
    echo "Building doc: $file"
    doxygen $file
done
