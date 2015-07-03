Integration test for protocols::helical_bundle::MakeBundle mover (which calls protocols::helical_bundle::MakeBundleHelix mover) when
the repeating unit in each helix is more than one residue.  This also tests the RosettaScripts connection.  The test creates a
collagen triple helix.  The structure is more or less meaningless, but it confirms the proper functioning of the mover. In this case,
the backbone generated is non-ideal (i.e. the mainchain bond angles and bond lengths are permitted to deviate a little bit from
ideality in order to generate a more perfect helix of helices).

This test also uses the PerturbBundle mover to perturb the generated helix, so it tests that, too.

Changes to the logfile or to the output pdb mean that the mover has changed somehow.  There should be no stochastic behaviour
associated with this mover, though numerical imprecision is a possibility.  Note that scores are expected to be extremely high,
since there is no energy-minimization or whatnot; this is just a crude backbone generation script.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).
