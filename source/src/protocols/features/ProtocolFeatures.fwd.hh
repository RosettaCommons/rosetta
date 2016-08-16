// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProtocolFeatures.fwd.hh
/// @brief  report Orbital geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ProtocolFeatures_fwd_hh
#define INCLUDED_protocols_features_ProtocolFeatures_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace features {

class ProtocolFeatures;
typedef utility::pointer::shared_ptr< ProtocolFeatures > ProtocolFeaturesOP;
typedef utility::pointer::shared_ptr< ProtocolFeatures const > ProtocolFeaturesCOP;

} // features namespace
} // protocols namespace

#endif //INCLUDED_protocols_features_ProtocolFeatures_fwd_hh
