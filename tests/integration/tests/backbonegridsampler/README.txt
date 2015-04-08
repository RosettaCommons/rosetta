Integration test for protocols::helical_bundle::BackboneGridSampler mover (which calls the
protocols::cyclic_peptide::PeptideStubMover mover).  This also tests the RosettaScripts connection.
The test creates a poly-ALA chain and samples its mainchain torsion space.  This structure is meaningless,
but it confirms the proper functioning of the mover.

Changes to the logfile or to the output pdb mean that the mover has changed somehow.  There should be
no stochastic behaviour associated with this mover, though numerical imprecision is a possibility.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).
