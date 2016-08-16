// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/MC_OneTorsion.cc
/// @brief Markov chain sampler for one torsion.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/MC_OneTorsion.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.MC_OneTorsion" );

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
MC_OneTorsion::MC_OneTorsion(
	core::id::TorsionID const & tor_id,
	core::Real const start_torsion
):
	MC_StepWiseSampler(),
	stored_angle_( start_torsion ),
	active_angle_( start_torsion ),
	angle_min_( -180 ),
	angle_max_( 180 ),
	stdev_( 10 ),
	torsion_id_( tor_id )
{}
///////////////////////////////////////////////////////////////////////////
void MC_OneTorsion::operator++() {
	if ( uniform_modeler() ) {
		active_angle_ = numeric::random::rg().uniform() * (angle_max_ - angle_min_) + angle_min_;
		regularize_angle( active_angle_ );
	} else {
		active_angle_ = numeric::random::rg().gaussian() * stdev_ + stored_angle_;
		if ( check_angle_in_range( active_angle_ ) ) {
			regularize_angle( active_angle_ );
		} else {
			active_angle_ = stored_angle_;
		}
	}
}
///////////////////////////////////////////////////////////////////////////
void MC_OneTorsion::apply( pose::Pose & pose ) {
	pose.set_torsion( torsion_id_, active_angle_ );
}
///////////////////////////////////////////////////////////////////////////
bool MC_OneTorsion::check_angle_in_range( Real angle ) const {
	if ( angle_max_ - angle_min_ >= 360 ) return true;
	while ( angle < angle_min_ ) angle += 360;
	while ( angle >= angle_max_ ) angle -= 360;
	return ( angle >= angle_min_ );
}
///////////////////////////////////////////////////////////////////////////
void MC_OneTorsion::regularize_angle( Real & angle ) {
	while ( angle < -180 ) angle += 360;
	while ( angle >= 180 ) angle -= 360;
}
///////////////////////////////////////////////////////////////////////////
} //sampler
} //stepwise
} //protocols

