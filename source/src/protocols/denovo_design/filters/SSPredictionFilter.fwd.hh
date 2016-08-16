// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/filters/SSPredictionFilter.fwd.hh
/// @brief header file for filter to determine agreement with psipred for secondary structure prediction
/// @details
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_filters_sspredictionfilter_fwd_hh
#define INCLUDED_protocols_denovo_design_filters_sspredictionfilter_fwd_hh

// Unit Headers

// Package headers

// Project headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace denovo_design {
namespace filters {

// Forward
class SSPredictionFilter;

// Types
typedef utility::pointer::shared_ptr< SSPredictionFilter > SSPredictionFilterOP;
typedef utility::pointer::shared_ptr< SSPredictionFilter const > SSPredictionFilterCOP;

}
}
}
#endif
