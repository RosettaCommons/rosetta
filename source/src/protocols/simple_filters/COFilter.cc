// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file COFilter.cc
/// @brief runs reject or accept filters on pose
/// @details
///	  Contains currently: COFilter
///
///
/// @author Robert Vernon

// Unit Headers
#include <protocols/simple_filters/COFilter.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/ContactOrderEnergy.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//// C++ headers
static thread_local basic::Tracer tr( "protocols.filters.COFilter" );

namespace protocols {
namespace simple_filters {

bool COFilter::apply( core::pose::Pose const & pose ) const{

	using core::Real;
	using core::Size;

	if ( !AbinitioBaseFilter::apply( pose ) ) return false;
	if ( sstype_ == "fail" ) return true;


	core::Real co_cutoff = 0.0;

	if ( sstype_ == "ab" ) {
		co_cutoff = ( 0.137f * pose.total_residue() ) + 3.25f;
	} else if ( sstype_ == "b" ) {
		co_cutoff = ( 0.145f * pose.total_residue() ) + 7.50f;
	} else {
		return true;
	}

	//Real co_score = 0.0;

  core::scoring::methods::ContactOrderEnergy co_energy;
  //core::scoring::EnergyMap emap;
  //core::scoring::ScoreFunction sfxn;

  //co_energy.finalize_total_energy( pose, sfxn, emap );
  Real co_score = co_energy.calculate_contact_order( pose );

	if ( co_score < co_cutoff ) {
		return false;
	}

	return true;
} // apply_filter

} // filters
} // protocols
