// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Phil Bradley


#include <protocols/rigid/rigid_body_moves.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>


//Utility Headers
#include <numeric/random/random.hh>


#include <basic/Tracer.hh>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;
namespace protocols {
namespace rigid {

using namespace core;


static numeric::random::RandomGenerator RG(6245211); // <- Magic number, do not change it!!!


int
gaussian_jump_move(
	pose::Pose & pose,
	kinematics::MoveMap const & mm,
	Real const translation_magnitude,
	Real const rotation_magnitude,
	int const dir // = 0 --> choose randomly
)
{

	utility::vector1< int > moving_jumps;
	for ( Size i=1, i_end = pose.num_jump(); i<= i_end; ++i ) {
		if ( mm.get_jump(i) ) moving_jumps.push_back( i );
	}

	if ( moving_jumps.empty() ) {
		T("protocols.rigid.rigid_body") << "[WARNING] no movable jumps!" << std::endl;
		return 0;
	}

	int const jump_number( RG.random_element( moving_jumps ) );
	gaussian_jump_move( pose, jump_number, translation_magnitude, rotation_magnitude, dir );
	return jump_number;
}

void
gaussian_jump_move(
	pose::Pose & pose,
	int const jump_number,
	Real const translation_magnitude,
	Real const rotation_magnitude,
	int dir // = 0 --> choose randomly
)
{
	kinematics::Jump jump( pose.jump( jump_number ) );
	if ( dir == 0 ) {
		dir = ( numeric::random::uniform() < 0.5 ? -1 : 1 );
	} else {
		runtime_assert( dir == 1 || dir == -1 );
	}
	jump.gaussian_move( dir, translation_magnitude, rotation_magnitude );
	pose.set_jump( jump_number, jump );
}



} // namespace rigid
} // namespace core
