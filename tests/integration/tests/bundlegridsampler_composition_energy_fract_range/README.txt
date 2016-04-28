Integration test for protocols::helical_bundle::BundleGridSampler mover (which calls the protocols::helical_bundle::MakeBundle mover), and of the aa_composition energy function (which is intended to guide packing to desired amino acid compositions).  This also tests the RosettaScripts connection.  This version tests a composition file that uses the FRACT_DELTA_START and FRACT_DELTA_END keywords.

The test creates a two-helix bundle and samples its parameter space, starting with helices too close together (which would normally give polyalaine sequences) and ending with helices too far apart.  This structure is meaningless, but it confirms the proper functioning of the mover.  In this case, the backbone generated is non-ideal (i.e. the mainchain bond angles and bond lengths are permitted to deviate a little bit from ideality in order to generate a more perfect helix of helices).

Changes to the logfile or to the output pdb mean that the mover has changed somehow.  There should be
no stochastic behaviour associated with this mover, though numerical imprecision is a possibility.

Note: when Alex Ford's changes to the packer that permit non-pairwise-decomposable score terms to be used in packing are implemented, this test should change, and should yield better sequences (with fewer alanines) in the cases of the too-close helices.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

