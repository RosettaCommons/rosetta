# pymol_cif description

This is a test to make sure that the CIF format that's output by PyMOL can be read in by Rosetta.

It's not 100% working, though, as the CIF library Rosetta uses will ignore the last table in the CIF.
(We get around this by adding a dummy table to the file, first.)
