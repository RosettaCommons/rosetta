Integration test for protocols::helical_bundle::MakeBundle mover (which calls protocols::helical_bundle::MakeBundleHelix mover).  This also tests the RosettaScripts connection.
The test creates a mixed alpha/beta structure consisting of an alpha-helical bundle surrounded by an antiparallel beta-barrel "fence".  The structure is more or less
meaningless, but it confirms the proper functioning of the mover.

Changes to the logfile or to the output pdb mean that the mover has changed somehow.  There should be no stochastic behaviour associated with this mover, though numerical
imprecision is a possibility.  Note that scores are expected to be extremely high, since there is no energy-minimization or whatnot; this is just a crude backbone generation
script.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

