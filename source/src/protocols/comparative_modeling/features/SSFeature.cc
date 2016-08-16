// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <protocols/comparative_modeling/features/ResidueFeature.hh>
#include <protocols/comparative_modeling/features/SSFeature.hh>
#include <protocols/comparative_modeling/features/SSFeature.fwd.hh>
#include <protocols/comparative_modeling/features/ResidueFeature.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {
namespace features {

SSFeature::SSFeature() :
	ss_type_( INVALID_SS )
{}

SSFeature::SSFeature( SSFeature const & other ) :
	ResidueFeature(),
	ss_type_( other.ss_type() )
{}

SSFeature::SSFeature( SSType ss ) :
	ss_type_( ss )
{}

utility::vector1< ResidueFeatureOP >
SSFeature::values_from_pose( core::pose::Pose & pose ) const {
	using core::Real;
	using std::string;
	using utility::vector1;

	vector1< SSFeatureOP > features;

	//protocols::jumping::assign_ss_dssp( pose );
	string const & secstruct( pose.secstruct() );
	for ( Size ii = 1; ii <= secstruct.size(); ++ii ) {
		features.push_back(
			protocols::comparative_modeling::features::SSFeatureOP( new SSFeature( char2ss_type(secstruct.at(ii)) ) )
		);
	}

	return features;
}

std::string SSFeature::type() const {
	return "ss";
}

ResidueFeatureOP SSFeature::clone() const {
	SSFeatureOP copy( new SSFeature( *this ) );
	return copy;
}

SSType SSFeature::ss_type() const {
	return ss_type_;
}

SSType SSFeature::char2ss_type( char const ss ) {
	if ( ss == 'H' ) {
		return H_SS;
	} else if ( ss == 'E' ) {
		return E_SS;
	} else if ( ss == 'L' ) {
		return L_SS;
	}

	return INVALID_SS;
}

} // features
} // comparative_modeling
} // protocols
