# Integration test "simple\_cycpep\_predict\_nterm\_isopeptide\_lariat"

##Author

Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

## Description

Integration test for the simple\_cycpep\_predict pilot app when used to close a macrocyclic
peptide cyclized by an isopeptide bond from an ASP or GLU side-chain to the N-terminus.  Failure
of this test means that something has changed about the way in which simple\_cycpep\_predict
handles isopeptide-cyclized peptides.  This version has a couple of "tail residues" which follow
the upper cyclization point.

## Expected output

Five PDB files containing isopeptide-cyclized peptides.  A logfile with RMSD to native.
