-----------------------
Notes on pdb_components
-----------------------
Rhiju Das, 2017

# How these files were prepared

1. Download components.cif from the PDB:

https://www.wwpdb.org/data/ccd

2. gunzip the file. Manually split the file into two pieces (A-C, D-Z), and gzip again; to get the file size under 25 Mb.

3. Update src/basic/options/options_rosetta.py PDB_components_file flag to point to the new files.


# TODO

A better strategy would be to keep the file uncompressed and let Git handle updates in a smart way. See discussion here:

https://github.com/RosettaCommons/main/pull/2693

One possibility would be to write a python script that ftps the current components.cif.gz and splits into:

 components.A.cif
 components.B.cif
  ...

And then have -PDB_components_file check through each one.

That could also save runtime if GlobalResidueTypeSet.cc only checks components.N.cif if it needs to find a name3 that starts with N.



