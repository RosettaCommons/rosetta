// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/McOneTorsion.cc
/// @brief Markov chain sampler for one torsion.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/McOneTorsion.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static basic::Tracer TR( "protocols.rotamer_sampler.McOneTorsion" );
static numeric::random::RandomGenerator RG( 2558493 );  // Magic number

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
McOneTorsion::McOneTorsion(
		core::id::TorsionID const & tor_id,
		core::Real const start_torsion
):
	McRotamer(),
	stored_angle_( start_torsion ),
	active_angle_( start_torsion ),
	angle_min_( -180 ),
	angle_max_( 180 ),
	stdev_( 10 ),
	torsion_id_( tor_id )
{}
///////////////////////////////////////////////////////////////////////////
void McOneTorsion::operator++() {
	if ( uniform_sampling() ) {
		active_angle_ = RG.uniform() * (angle_max_ - angle_min_) + angle_min_;
		regularize_angle( active_angle_ );
	} else {
		active_angle_ = RG.gaussian() * stdev_ + stored_angle_;
		if ( check_angle_in_range( active_angle_ ) ) {
			regularize_angle( active_angle_ );
		} else {
			active_angle_ = stored_angle_;
		}
	}
}
///////////////////////////////////////////////////////////////////////////
void McOneTorsion::apply( pose::Pose & pose ) {
	pose.set_torsion( torsion_id_, active_angle_ );
}
///////////////////////////////////////////////////////////////////////////
bool McOneTorsion::check_angle_in_range( Real angle ) const {
	if ( angle_max_ - angle_min_ >= 360 ) return true;
	while ( angle < angle_min_ ) angle += 360;
	while ( angle >= angle_max_ ) angle -= 360;
	return ( angle >= angle_min_ );
}
///////////////////////////////////////////////////////////////////////////
void McOneTorsion::regularize_angle( Real & angle ) {
	while ( angle < -180 ) angle += 360;
	while ( angle >= 180 ) angle -= 360;
}
///////////////////////////////////////////////////////////////////////////
} //rotamer_sampler
} //protocols

