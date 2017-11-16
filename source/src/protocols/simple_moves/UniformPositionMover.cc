// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/simple_moves/UniformPositionMover.hh
///
/// @brief      Apply a uniform (deterministic move) to a given position
/// @details Generic movers for applying a uniform rotation or translation move
///    along a jump to a given partner in the pose.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/10/14)

// Unit Headers
#include <protocols/simple_moves/UniformPositionMover.hh>

// Project Headers
#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>

#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

static basic::Tracer TR( "protocols.simple_moves.UniformPositionMover" );

namespace protocols {
namespace simple_moves {

using namespace core;

// Position Rotation Mover /////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Uniform Position Mover
UniformRotationMover::UniformRotationMover() :
	Mover(),
	alpha_( 0.0 ),
	axis_( 0.0, 0.0, 1.0 ),
	rb_jump_( 1 )
{}

/// @brief Custom Constructor
/// @details Specify a new normal axis
UniformRotationMover::UniformRotationMover(
	Real alpha,
	Vector axis,
	core::SSize rb_jump
) :
	Mover(),
	alpha_( alpha ),
	axis_( axis ),
	rb_jump_( rb_jump )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
UniformRotationMover::UniformRotationMover( UniformRotationMover const & ) = default;

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
UniformRotationMover &
UniformRotationMover::operator=( UniformRotationMover const & src )
{

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new UniformRotationMover( *this ) );

}

/// @brief Destructor
UniformRotationMover::~UniformRotationMover() = default;

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
UniformRotationMover::get_name() const {
	return "UniformRotationMover";
}

/// @brief Apply Rotation move
/// @brief Rotate position over jump given an axis & angle
void
UniformRotationMover::apply( Pose & pose ) {

	using namespace numeric;
	using namespace core::kinematics;

	// Grab the upstream & downstream stubs from the pose
	Stub downstream_stub( pose.conformation().downstream_jump_stub( rb_jump_ ) );
	Stub upstream_stub( pose.conformation().upstream_jump_stub( rb_jump_ ) );

	// Compute Rotation (delta_rot)
	xyzMatrix< core::Real > delta_rot = rotation_matrix( axis_, alpha_ );
	downstream_stub.M = delta_rot * downstream_stub.M;

	Jump flexible_jump = Jump( RT( downstream_stub, upstream_stub ) );

	// Set jump in the pose
	pose.set_jump( rb_jump_, flexible_jump );

}


// Position Translation Mover //////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Translation Mover
UniformTranslationMover::UniformTranslationMover() :
	Mover(),
	new_position_( 0.0, 0.0, 0.0 ),
	rb_jump_( 1 )
{}

/// @brief Custom Constructor
/// @details Specify a new point to move to
UniformTranslationMover::UniformTranslationMover(
	Vector new_position,
	SSize rb_jump
) :
	Mover(),
	new_position_( new_position ),
	rb_jump_( rb_jump )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
UniformTranslationMover::UniformTranslationMover( UniformTranslationMover const & ) = default;

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
UniformTranslationMover &
UniformTranslationMover::operator=( UniformTranslationMover const & src )
{

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new UniformTranslationMover( *this ) );

}

/// @brief Destructor
UniformTranslationMover::~UniformTranslationMover() = default;

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
UniformTranslationMover::get_name() const {
	return "UniformTranslationMover";
}

/// @brief Apply Translation
void
UniformTranslationMover::apply( Pose & pose ) {

	using namespace numeric;
	using namespace core::kinematics;

	// Grab the upstream and downstream stubs in the foldtree
	Stub downstream_stub( pose.conformation().downstream_jump_stub( rb_jump_ ) );
	Stub upstream_stub( pose.conformation().upstream_jump_stub( rb_jump_ ) );

	downstream_stub.v = new_position_;

	// Create a jump between the upstream and downsteam jump
	Jump flexible_jump = Jump( RT( upstream_stub, downstream_stub ) );

	// Set jump in the pose
	pose.set_jump( rb_jump_, flexible_jump );

}


} // simple_moves
} // protocols
