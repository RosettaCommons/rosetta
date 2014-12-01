Integration test for protocols::helical_bundle::PerturbBundle mover (which calls protocols::helical_bundle::PerturbBundleHelix mover).  This also tests the RosettaScripts connection.
The test creates a mixed alpha/beta structure consisting of an alpha-helical bundle surrounded by an antiparallel beta-barrel "fence".  The structure is more or less
meaningless, but it confirms the proper functioning of the mover.  The bundle is then perturbed randomly, subject to certain user settings, by the PerturbBundle mover.

In this case, the backbone generated is non-ideal (i.e. the mainchain bond angles and bond lengths are permitted to deviate a little bit from ideality in order to generate a more perfect helix of helices).

Changes to the logfile or to the output pdb mean that the mover has changed somehow.  Note that the MakeBundle mover is deterministic, while the PerturbBundle mover is stocahstic.
Note also that scores are expected to be extremely high, since there is no energy-minimization or whatnot; this is just a crude backbone generation script.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

