# Integration test "simple\_cycpep\_predict\_cterm\_isopeptide\_lariat\_tailless"

##Author

Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

## Description

Integration test for the simple\_cycpep\_predict pilot app when used to close a macrocyclic
peptide cyclized by an isopeptide bond from a nitrogenous side-chain to the C-terminus.  Failure
of this test means that something has changed about the way in which simple\_cycpep\_predict
handles isopeptide-cyclized peptides.  This version has a couple of "tail residues" which precede
the lower cyclization point.

This version is for a cyclic peptide with no tail.

## Expected output

Five PDB files containing isopeptide-cyclized peptides.
