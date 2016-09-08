// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RGFilter.cc
/// @brief runs reject or accept filters on pose
/// @details
///   Contains currently: RGFilter
///
///
/// @author Robert Vernon

// Unit Headers
#include <protocols/simple_filters/RGFilter.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//// C++ headers
static THREAD_LOCAL basic::Tracer tr( "protocols.simple_filters.RGFilter" );

namespace protocols {
namespace simple_filters {

bool RGFilter::apply( core::pose::Pose const & pose ) const {

	using namespace core::scoring;
	if ( !AbinitioBaseFilter::apply( pose ) ) return false;
	if ( sstype_ == "fail" ) return true;
	//car sheet_filter disable

	//car rg filter disable
	if ( max_helix_length_ > 25 || max_helix_fraction_ > 0.22 ) return true;

	//ScoreFunctionOP scorefxn( new ScoreFunction );
	//scorefxn->set_weight( rg, 1.0 );
	core::scoring::methods::RG_Energy_Fast rg_energy;

	core::Real rg_threshold
		= ( 3.0 * std::pow( static_cast< float >( pose.size() ), ( 1.0f / 3 ) ) ) + 2;
	//core::Real rg_score = (*scorefxn)( pose );

	core::Real rg_score = rg_energy.calculate_rg_score( pose );


	if ( rg_score > rg_threshold ) {
		return false;
	}

	return true;
}

} // filters
} // protocols
