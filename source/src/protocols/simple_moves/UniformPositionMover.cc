// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/simple_moves/UniformPositionMover.hh
///
/// @brief      Apply a uniform (deterministic move) to a given position
/// @details	Generic movers for applying a uniform rotation or translation move
///				along a jump to a given partner in the pose.
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
UniformPositionRotationMover::UniformPositionRotationMover() :
	Mover(),
	alpha_( 0.0 ),
	axis_( 0, 0, 1 ),
	rb_jump_( 1 )
{}

/// @brief Custom Constructor
/// @details Specify a new normal axis
UniformPositionRotationMover::UniformPositionRotationMover(
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
UniformPositionRotationMover::UniformPositionRotationMover( UniformPositionRotationMover const & src ) :
	Mover( src ),
	alpha_( src.alpha_ ),
	axis_( src.axis_ ),
	rb_jump_( src.rb_jump_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
UniformPositionRotationMover &
UniformPositionRotationMover::operator=( UniformPositionRotationMover const & src )
{
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new UniformPositionRotationMover( *this ) );
	
}

/// @brief Destructor
UniformPositionRotationMover::~UniformPositionRotationMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
UniformPositionRotationMover::get_name() const {
	return "UniformPositionRotationMover";
}

/// @brief Apply Rotation move
/// @brief Rotate position over jump given an axis & angle
void
UniformPositionRotationMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	
	// Grab the upstream & downstream stubs from the pose
	Stub membrane_stub( pose.conformation().downstream_jump_stub( rb_jump_ ) );
	Stub anchor_stub( pose.conformation().upstream_jump_stub( rb_jump_ ) );
	
	// Compute Rotation (detla_rot)
	xyzMatrix< core::Real > delta_rot = rotation_matrix( axis_, alpha_ );
	membrane_stub.M = delta_rot * membrane_stub.M;
	
	Jump flexible_jump = Jump( RT( membrane_stub, anchor_stub ) );
	
	// Set jump in the pose
	pose.set_jump( rb_jump_, flexible_jump );
	
}


// Position Translation Mover /////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
UniformPositionTranslationMover::UniformPositionTranslationMover() :
	Mover(),
	new_position_( 0, 0, 0 ),
	rb_jump_( 1 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane center to move to
UniformPositionTranslationMover::UniformPositionTranslationMover(
	Vector new_position,
	SSize rb_jump
	) :
	Mover(),
	new_position_( new_position ),
	rb_jump_( rb_jump )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
UniformPositionTranslationMover::UniformPositionTranslationMover( UniformPositionTranslationMover const & src ) :
	Mover( src ),
	new_position_( src.new_position_ ),
	rb_jump_( src.rb_jump_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
UniformPositionTranslationMover &
UniformPositionTranslationMover::operator=( UniformPositionTranslationMover const & src )
{
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new UniformPositionTranslationMover( *this ) );
	
}

/// @brief Destructor
UniformPositionTranslationMover::~UniformPositionTranslationMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
UniformPositionTranslationMover::get_name() const {
	return "UniformPositionTranslationMover";
}

/// @brief Apply Rotation/Translation to Membrane
/// @brief Translate the membrane position in this pose
/// to the new center position, and rotate to new normal
void
UniformPositionTranslationMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	
	// Grab the stub of the membrane residue & its anchoring jump
	Stub downstream_stub( pose.conformation().downstream_jump_stub( rb_jump_ ) );
	Stub upstream_stub( pose.conformation().upstream_jump_stub( rb_jump_ ) );
	
	downstream_stub.v = new_position_;
	
	// Create a jump between the membrane and its anchoring point
	Jump flexible_jump = Jump( RT( upstream_stub, downstream_stub ) );
	
	// Set jump in the pose
	pose.set_jump( rb_jump_, flexible_jump );
	
}


} // simple_moves
} // protocols
