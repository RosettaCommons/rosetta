// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ResidueTotalScoresFeatures.fwd.hh
/// @brief  report the total per residue score to a features database
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ResidueTotalScoresFeatures_fwd_hh
#define INCLUDED_protocols_features_ResidueTotalScoresFeatures_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace features {

class ResidueTotalScoresFeatures;
typedef utility::pointer::shared_ptr< ResidueTotalScoresFeatures > ResidueTotalScoresFeaturesOP;
typedef utility::pointer::shared_ptr< ResidueTotalScoresFeatures const > ResidueTotalScoresFeaturesCOP;

}//features
}//protocols

#endif //INCLUDED_protocols_features_ResidueTotalScoresFeatures_fwd_hh
