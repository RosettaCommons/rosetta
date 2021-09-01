# Integration test "helical\_bundle\_predict\_psipred\_with\_helix\_globals"

##Author

Vikram K. Mulligan, Ph.D., Center for Computational Biology, Flatiron Institute (vmulligan@flatironinstitute.org).

## Description

Integration test for the helical\_bundle\_predict application, with a PsiPred prediction file plus a helix definition
file that sets global sampling parameters.  This tests prediction of PDB file 1b67.  The test is not expected to produce
a particularly accurate prediction with an nstruct of 5.  The main thing is to ensure that this runs. 

## Expected output

Five PDB files containing vaguely helical bundle-like polypeptides.
