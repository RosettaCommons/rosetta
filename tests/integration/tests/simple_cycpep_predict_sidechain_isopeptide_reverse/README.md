# Integration test "simple\_cycpep\_predict\_sidechain\_isopeptide\_reverse"

##Author

Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

## Description

Integration test for the simple\_cycpep\_predict pilot app when used to close a peptide
macrocycle that's cyclized through a side-chain isopeptide bond.  Failure of this test
means that something has changed about the way in which simple\_cycpep\_predict handles
isopeptide-cyclized peptides.  This version has a few "tail residues" which precede the
lower cyclization point and follow the upper.

This version has the carboxyl and amide residues reversed compared to the
simple\_cycpep\_predict\_sidechain\_isopeptide integration test.

## Expected output

Five PDB files containing sidechain isopeptide-cyclized peptides.  There should be some
tail residues at each end.
