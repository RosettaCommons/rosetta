// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/constraints/ResidueFeatureFactory.hh
/// @brief Factory for creating various types of constraints.
/// @author James Thompson <tex@u.washington.edu>

#ifndef INCLUDED_protocols_comparative_modeling_features_ResidueFeatureFactory_hh
#define INCLUDED_protocols_comparative_modeling_features_ResidueFeatureFactory_hh

// Package headers
#include <protocols/comparative_modeling/features/ResidueFeature.fwd.hh>

// Project headers

// C++ Headers
#include <map>

namespace protocols {
namespace comparative_modeling {
namespace features {

class ResidueFeatureFactory {
public:
	ResidueFeatureFactory(void);

 	/// @brief adds a ResidueFeatureOP
	void add_type( ResidueFeatureOP feature_ );
	ResidueFeatureOP get_residue_feature( std::string const & type ) const;

	typedef std::map< std::string, ResidueFeatureOP > ResidueFeatureTypes;
	ResidueFeatureTypes feature_types_;
};

} // features
} // sequence
} // core

#endif
