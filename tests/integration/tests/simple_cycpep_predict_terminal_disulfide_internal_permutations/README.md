# Integration test "simple\_cycpep\_predict\_terminal\_internal\_permutations

##Author

Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

## Description

Integration test for the simple\_cycpep\_predict pilot app when used to close a disulfize-
cyclized peptide.  Failure of this test means that something has changed about the way
in which simple\_cycpep\_predict handles disulfide-cyclized peptides.  This version tries
all permutations of internal disulfides, too.

## Expected output

Five PDB files containing disulfide-cyclized poly-gly peptides.  The terminal disulfide
should be between a DCYS and an LCYS.  There should be two tail residues at the N-terminus,
and none at the C-terminus.  Two internal disulfides should also be present.
