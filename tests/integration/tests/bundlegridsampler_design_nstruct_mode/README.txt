This test is intended to test the integration of the BundleGridSampler mover with the FastDesign mover and the LayerDesign TaskOperation in nstruct mode (splitting Crick parameter samples over jobs specified with the -nstruct command).  More generally, it tests the ability to pass other movers to the BundleGridSampler mover in RosettaScripts.  The test has now been updated so that it also tests the ability to pass filters to the BundleGridSampler mover in RosettaScripts.

Note that although the BundleGridSampler is not stochastic, the FastDesign mover is, so that might
be an issue.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

