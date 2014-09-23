// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/features/ResidueFeatureFactory.cc
/// @brief Factory for creating various types of ResidueFeatures.
/// @author James Thompson

// Unit headers
#include <protocols/comparative_modeling/features/ResidueFeature.hh>
#include <protocols/comparative_modeling/features/ResidueFeature.fwd.hh>
#include <protocols/comparative_modeling/features/ResidueFeatureFactory.hh>

// Package headers
#include <protocols/comparative_modeling/features/SSFeature.hh>
#include <protocols/comparative_modeling/features/TorsionFeature.hh>

#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {
namespace features {

void
ResidueFeatureFactory::add_type( ResidueFeatureOP new_feature ) {
	std::string const type_name( new_feature->type() );
	feature_types_[ type_name ] = new_feature;
}

ResidueFeatureOP
ResidueFeatureFactory::get_residue_feature(
	std::string const & type
) const {
	ResidueFeatureTypes::const_iterator iter = feature_types_.find( type );
	if ( iter != feature_types_.end() ) {
		return iter->second->clone();
	} else {
		utility_exit_with_message(
			"ResidueFeatureFactory: unknown ResidueFeature type: " + type
		);
		return NULL;
	}
}

ResidueFeatureFactory::ResidueFeatureFactory(void) {
	// initialization of ResidueFeatures which this factory knows how to
	// instantiate
	add_type( ResidueFeatureOP( new SSFeature() ) );
	add_type( ResidueFeatureOP( new TorsionFeature() ) );
}

} // namespace features
} // namespace comparative_modeling
} // namespace core
