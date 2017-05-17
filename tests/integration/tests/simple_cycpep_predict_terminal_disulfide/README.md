# Integration test "simple\_cycpep\_predict\_terminal\_disulfide"

##Author

Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

## Description

Integration test for the simple\_cycpep\_predict pilot app when used to close a disulfize-
cyclized peptide.  Failure of this test means that something has changed about the way
in which simple\_cycpep\_predict handles disulfide-cyclized peptides.

## Expected output

Five PDB files containing disulfide-cyclized poly-gly peptides.  The terminal disulfide
should be between a DCYS and an LCYS.
