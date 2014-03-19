#!/bin/bash

echo 'Building documentation for each src/doxygen.*'

./update_options.sh

mkdir html
cp src/index.html html/

rm -rf html/core+protocols html/all_else

for file in src/Doxyfile.* ; do
    echo "Building doc: $file"
    doxygen $file
done
