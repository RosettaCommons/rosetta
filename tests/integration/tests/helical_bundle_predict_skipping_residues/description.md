# Integration test "helical\_bundle\_predict\_skipping\_residues"

##Author

Vikram K. Mulligan, Ph.D., Center for Computational Biology, Flatiron Institute (vmulligan@flatironinstitute.org).

## Description

Integration test for the helical\_bundle\_predict application.  This attempts to predict the structure of Chris Bahl's miniprotein
with PDB ID 2ND2.  Note that it is not expected that the structure will be correct in the short time of an integration test.  The
main thing is to ensure that the app runs.

This version tests flags that allow poses of unequal length to be compared.  Two additional glycines that are not present in the
native structure have been inserted in the first loop, and these are ignored in the alignment, as are the first two and last two
residues of the pose.

## Expected output

Five PDB files containing vaguely helical bundle-like polypeptides.
