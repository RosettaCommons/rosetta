Integration test for protocols::helical_bundle::MakeBundle mover (which calls protocols::helical_bundle::MakeBundleHelix mover).  This also tests the RosettaScripts connection.
The test creates a mixed alpha/beta structure consisting of an alpha-helical bundle surrounded by an antiparallel beta-barrel "fence".  The structure is more or less
meaningless, but it confirms the proper functioning of the mover.  In this case, the backbone generated is non-ideal (i.e. the mainchain bond angles and bond lengths are permitted to deviate a little bit from ideality in order to generate a more perfect helix of helices).

Changes to the logfile or to the output pdb mean that the mover has changed somehow.  There should be no stochastic behaviour associated with this mover, though numerical
imprecision is a possibility.  Note that scores are expected to be extremely high, since there is no energy-minimization or whatnot; this is just a crude backbone generation
script.

Test 2 has been added to test the mover with input in degrees instead of radians.  The generated structure should be identical (within numerical error) to the one created using radians.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

