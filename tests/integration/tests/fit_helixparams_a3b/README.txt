Integration test for protocols::helical_bundle::FitSimpleHelix mover and the apps/pilot/vmullig/fit_helixparams.cc app.
 
The test creates a standard alpha-helix from L-alanine, and fits the minor helix Crick parameters to it, writing out the fitted parameters.

Changes to the logfile or to the output pdb / .crick_params file mean that the fitter has changed.  There should be no stochastic behaviour to the app or the mover, though numerical imprecision during nonlinear least-squares fitting (which uses the minimizer) may be an issue.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

