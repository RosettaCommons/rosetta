Integration test for protocols::simple_moves::SetTorsion mover.  This tests the setting, perturbing, randomizing, and randomizing by Rama functionality of the mover.  Failure of this test means that the SetTorsion class has changed, or that there's some issue with the stability of the random number generators.

This test also uses the PeptideStubMover.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).
