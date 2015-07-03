Integration tests for the BundleGridSampler mover for helices in which the repeating unit is more than one residue (in this case, the collagen triple-helix).  This generates a triple helix of collagen triple helices.  The structure is pretty meaningless, but it tests that the mover works properly with RosettaScripts and everything else.  Chages to this test could mean:

--RosettaScripts has changed somehow.
--The BundleGridSampler has changed.
--The MakeBundle or MakeBundleHelix movers have changed.
--Parameter storage in poses or BundleParameters have changed.
--The collagen.crick_params file in the database has changed.

The test should have no stochastic component to its behaviour.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).
