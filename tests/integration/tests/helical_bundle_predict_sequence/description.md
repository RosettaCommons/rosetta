# Integration test "helical\_bundle\_predict"

##Author

Vikram K. Mulligan, Ph.D., Center for Computational Biology, Flatiron Institute (vmulligan@flatironinstitute.org).

## Description

Integration test for the helical\_bundle\_predict application.  This attempts to predict the structure of Chris Bahl's miniprotein
with PDB ID 2ND2.  Note that it is not expected that the structure will be correct in the short time of an integration test.  The
main thing is to ensure that the app runs.

This version checks full-sequence input, and includes a GLY to DSER mutation.

## Expected output

Five PDB files containing vaguely helical bundle-like polypeptides.
