-----------------------
Notes on pdb_components
-----------------------
Rhiju Das, 2017

# How these files were prepared

1. Download components.cif from the PDB:

https://www.wwpdb.org/data/ccd

2. gunzip the file. Manually split the file into two pieces (A-C, D-Z), and gzip again; to get the file size under 25 Mb.

3. Update src/basic/options/options_rosetta.py PDB_components_file flag to point to the new files.

Andy Watkins, 2018

Rhiju's proposed 'better strategy' -- divide the components by initial letter and let git version them -- have been updated. If you want to update the components, run update_components.sh in this directory, which only requires curl and python.

