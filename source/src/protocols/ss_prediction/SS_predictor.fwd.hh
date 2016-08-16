// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ss_prediction/SS_predictor.fwd.hh
/// @brief  SS_predictor forward header
/// @author TJ Brunette


#ifndef INCLUDED_protocols_ss_prediction_SS_predictor_fwd_hh
#define INCLUDED_protocols_ss_prediction_SS_predictor_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace ss_prediction {

// Forward
class SS_predictor;

// Types
typedef  utility::pointer::shared_ptr< SS_predictor >  SS_predictorOP;
typedef  utility::pointer::shared_ptr< SS_predictor const >  SS_predictorCOP;

typedef  utility::pointer::weak_ptr< SS_predictor >  SS_predictorAP;
typedef  utility::pointer::weak_ptr< SS_predictor const >  SS_predictorCAP;


} // namespace ss_prediction
} // namespace protocols

#endif
