-----------------------
Notes on pdb_components
-----------------------

To update:
----------

Simply run update_components.sh. This will automatically download the current components.cif file from the wwPDB (using curl) and split it out into the respective alphabetized sub lists.

You should be able to simply commit and push the changes.

Overrides
---------

Some of the component definitions from the wwPDB are junk. (The crystallographers which submitted them did only a perfunctory job of it.)
Others are okay for the purposes of the PDB, but cause issues with the best-effort ResidueType generation system Rosetta has.
To account for this, we have an overrides/ directory. 

Simply put the fixed definition (one file per ligand) into the overrides/ directory, and then put the filename/path to the file in the overrides.txt file. 
