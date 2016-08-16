// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/BetaTurnDetection.fwd.hh
/// @brief  determine if there is a beta turn at particular postion
/// @author Brian Weitzner

#ifndef INCLUDED_protocols_features_BetaTurnDetection_FWD_HH
#define INCLUDED_protocols_features_BetaTurnDetection_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace features {

class BetaTurnDetection;
typedef utility::pointer::shared_ptr< BetaTurnDetection > BetaTurnDetectionOP;
typedef utility::pointer::shared_ptr< BetaTurnDetection const > BetaTurnDetectionCOP;

}//features
}//protocols

#endif //INCLUDED_protocols_features_BetaTurnDetection_FWD_HH
