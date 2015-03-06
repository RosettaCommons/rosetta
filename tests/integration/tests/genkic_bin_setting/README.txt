Integration test for protocols::generalized_kinematic_closure::GeneralizedKIC mover using core::scoring::bin_transitions::BinTransitionCalculator class to specify bins before kinematic closure (the set_backbone_bin GeneralizedKICperturber).  This test also tests the backbone_bin GeneralizedKICfilter.  Failure of this test means that the BinTransitionCalculator class has changed, that bin transition database files (located in database/protocol_data/generalizedKIC/bin_params/) have changed, or that GeneralizedKIC, the GeneralizedKICfilter, or the GeneralizedKICperturber have changed.

Author: Vikram K. Mulligan, Ph.D., Baker laboratory, University of Washington (vmullig@uw.edu).

