// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/ss_prediction/SS_predictorFromSilents.fwd.hh
/// @brief  SS_predictorFromSilents 
/// @author TJ Brunette


#ifndef INCLUDED_protocols_ss_prediction_SS_predictorFromSilents_fwd_hh
#define INCLUDED_protocols_ss_prediction_SS_predictorFromSilents_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace ss_prediction {

// Forward
class SS_predictorFromSilents;

// Types
typedef  utility::pointer::shared_ptr< SS_predictorFromSilents >  SS_predictorFromSilentsOP;
typedef  utility::pointer::shared_ptr< SS_predictorFromSilents const >  SS_predictorFromSilentsCOP;

typedef  utility::pointer::weak_ptr< SS_predictorFromSilents >  SS_predictorFromSilentsAP;
typedef  utility::pointer::weak_ptr< SS_predictorFromSilents const >  SS_predictorFromSilentsCAP;


} // namespace ss_prediction
} // namespace protocols

#endif
