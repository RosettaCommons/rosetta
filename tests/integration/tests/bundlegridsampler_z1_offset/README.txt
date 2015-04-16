Integration test for protocols::helical_bundle::BundleGridSampler mover (which calls the
protocols::helical_bundle::MakeBundle mover), testing the z1_offset sampling.  This also
tests the RosettaScripts connection. The test creates a four-helix bundle and samples its
parameter space.  This structure is meaningless, but it confirms the proper functioning
of the mover.  In this case, the backbone generated is non-ideal (i.e. the mainchain bond
angles and bond lengths are permitted to deviate a little bit from ideality
in order to generate a more perfect helix of helices).

Changes to the logfile or to the output pdb mean that the mover has changed somehow.  There should be
no stochastic behaviour associated with this mover, though numerical imprecision is a possibility.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

