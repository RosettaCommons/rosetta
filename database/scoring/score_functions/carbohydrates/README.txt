This directory contains lists of Gaussian parameters for the multiple terms of
the various CarboHydrate Intrinsic (CHI) energy functions, as defined in:
A.K. Nivedha et al. J. Comput. Chem. 2014, 35, 526-39
and
A.K. Nivedha et al. JCTC 2016, 12, 892-901

The parameters for energy functions in this directory are linked to a
LinkageType enum.  To add a new energy function, one must A) add a new
LinkageType to LinkageType.hh; B) add the filename of the parameters to the map
of enums to filenames in core::scoring::carbohydratesCHIEnergyFunction::init();
and C) add code to core::scoring::methods::carbohydrates::SugarBackboneEnergy.cc
to properly assign the new LinkageType.
