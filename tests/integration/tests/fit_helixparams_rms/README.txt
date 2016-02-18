Integration test for protocols::helical_bundle::FitSimpleHelix mover and the apps/pilot/vmullig/fit_helixparams.cc app.
 
The test creates am alpha-sheet from alternating L-alanine/D-alanine, and fits the minor helix Crick parameters to it, writing out the fitted parameters.  It tests that the RMS code is working properly -- I fixed this recently, and the code is now being refactored, so I want to make sure it stays fixed.

Changes to the logfile or to the output pdb / .crick_params file mean that the fitter has changed or that the rms code has changed.  There should be no stochastic behaviour to the app or the mover, though numerical imprecision during nonlinear least-squares fitting (which uses the minimizer) may be an issue.

The output PDB file should have copper atoms overlaid perfectly on all of the mainchain heavyatoms (except the carbonyl oxygens).  The output geometry should not be distorted.  INSPECT THIS OUTPUT MANUALLY TO CONIFRM THIS.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

