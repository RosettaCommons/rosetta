// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMLJEnergyTable.cc
/// @brief  Molecular mechanics lj energy table class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/mm/MMLJEnergyTable.hh>
#include <core/scoring/mm/MMLJScore.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.hh>


// Utility header
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <map>

#include <core/scoring/mm/MMLJLibrary.hh>
#include <utility/vector1.hh>
#include <cmath>


namespace core {
namespace scoring {
namespace mm {

/// @details Auto-generated virtual destructor
MMLJEnergyTable::~MMLJEnergyTable() {}

static THREAD_LOCAL basic::Tracer TR( "core.scoring.mm.MMLJEnergyTable" );

MMLJEnergyTable::MMLJEnergyTable() :
	max_dist_(49.0),
	linear_switch_point(0.6) // 60% of distance at minimum
{}

void
MMLJEnergyTable::score( Size atom1, Size atom2, Size path_distance, Real squared_distance, Real & rep, Real & atr ) const
{
	// init values
	rep = atr = 0;

	if ( squared_distance > max_dist_ ) return;
	if ( path_distance != 3 ) path_distance = 4;

	// subtract max dist slope and value
	// so when we hit that value, it's 0
	Real max_ener  = mm_lj_score_.score( atom1, atom2, path_distance, max_dist_ );
	Real max_slope = mm_lj_score_.deriv_score( atom1, atom2, path_distance, max_dist_ );
	Real max_int   = max_ener - max_slope * max_dist_;

	// get the distance and energy for when the function is a minimum
	Real min_ener_dist = mm_lj_score_.min_dist( atom1, atom2, path_distance );
	Real min_ener      = mm_lj_score_.score( atom1, atom2, path_distance, min_ener_dist*min_ener_dist );

	// get values for linear switching at short distances
	Real switch_dist         = linear_switch_point * min_ener_dist;
	Real switch_dist_squared = switch_dist * switch_dist;
	Real switch_slope = mm_lj_score_.deriv_score( atom1, atom2, path_distance, switch_dist_squared );
	Real switch_ener  = mm_lj_score_.score(       atom1, atom2, path_distance, switch_dist_squared );
	Real switch_intercept = switch_ener - switch_slope * switch_dist_squared;

	Real temp_score = 0;
	if ( squared_distance <= switch_dist_squared ) { // in switch region
		temp_score = switch_slope * squared_distance + switch_intercept;
	} else {
		temp_score = mm_lj_score_.score( atom1, atom2, path_distance, squared_distance );
	}

	if ( squared_distance < min_ener_dist * min_ener_dist ) { // repulsive
		rep = temp_score - min_ener;
		atr = min_ener - (max_int + max_slope * squared_distance );//max_ener;
	} else { // attractive
		rep = 0;
		atr = temp_score - (max_int + max_slope * squared_distance );//max_ener;
	}
}

void
MMLJEnergyTable::deriv_score( Size atom1, Size atom2, Size path_distance, Real squared_distance, Real & drep, Real & datr ) const
{
	// init values
	drep = datr = 0;
	if ( squared_distance > max_dist_ ) return;

	// Split out behavior for path distance 3 only
	if ( path_distance != 3 ) path_distance = 4;

	Real max_slope = mm_lj_score_.deriv_score( atom1, atom2, path_distance, max_dist_ );

	// get the distance and energy for when the function is a minimum
	Real min_ener_dist = mm_lj_score_.min_dist( atom1, atom2, path_distance );

	// get values for linear switching at short distances
	Real switch_dist = linear_switch_point * min_ener_dist;
	Real switch_dist_squared = switch_dist * switch_dist;
	Real switch_slope = mm_lj_score_.deriv_score( atom1, atom2, path_distance, switch_dist_squared );

	Real temp_deriv = 0;
	if ( squared_distance <= switch_dist_squared ) { // in switch region
		temp_deriv = switch_slope;
	} else {
		temp_deriv = mm_lj_score_.deriv_score( atom1, atom2, path_distance, squared_distance );
	}

	if ( squared_distance < min_ener_dist * min_ener_dist ) { // repulsive
		drep = temp_deriv;
		datr = 0 - max_slope;
	} else { // attractive
		drep = 0;
		datr = temp_deriv - max_slope;
	}
}

} // namespace mm
} // namespace scoring
} // namespace core
