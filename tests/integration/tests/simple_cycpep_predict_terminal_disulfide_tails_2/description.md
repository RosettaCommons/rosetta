# Integration test "simple\_cycpep\_predict\_terminal\_disulfide\_tails\_2"

##Author

Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org).

## Description

Integration test for the simple\_cycpep\_predict pilot app when used to close a disulfize-
cyclized peptide.  Failure of this test means that something has changed about the way
in which simple\_cycpep\_predict handles disulfide-cyclized peptides.  This version has a 
couple of "tail residues" which precede the lower cysteine and follow the upper.  It is a
second test case that the first test was not catching.

## Expected output

Five PDB files containing disulfide-cyclized poly-gly peptides.  The terminal disulfide
should be between two LCYS residues.  There should be one tail residues at the N-terminus
and two at the C-terminus.
