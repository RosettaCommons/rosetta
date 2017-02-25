This tests the RamaMutationSelector, which selects positions based on the RamaPrepro energy of a residue at that position IF
the position were mutated to a particular residue type.  This is useful for picking positions to mutate to a conformationally-constrained
type, like proline or AIB.  Changes to this test mean the selector has changed somehow, or possibly the MutateResidue mover, which
the test also calls.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

