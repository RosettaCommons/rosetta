// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/RandomMembranePositionMover.hh
///
/// @brief      Random membrane position mover
/// @details Make random perturbations in the position of the membrane
///    as define by some center position and normal coordinate
///    (also of fixed thickness). Rotation mover rotates normal
///    from a random angle, Translation mover picks a randomly
///    translated position from the current position.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/10/14)

#ifndef INCLUDED_protocols_membrane_RandomMembranePositionMover_hh
#define INCLUDED_protocols_membrane_RandomMembranePositionMover_hh

// Unit Headers
#include <protocols/membrane/RandomMembranePositionMover.fwd.hh>

// Projact Headers
#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace membrane {


/// @brief Random Position Rotation Move
/// @details Rotate the orientation of the membrane position to a new
/// normal position. Angle chosen randomnly - moved from current axis
class RandomPositionRotationMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Construct a Default Membrane Position Mover
	RandomPositionRotationMover();

	/// @brief Custom Constructor
	/// @details Specify an order of magnitude to rotate
	RandomPositionRotationMover(
		core::Real rot_mag,
		core::SSize rb_jump
	);

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	RandomPositionRotationMover( RandomPositionRotationMover const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	RandomPositionRotationMover &
	operator=( RandomPositionRotationMover const & src );

	/// @brief Destructor
	~RandomPositionRotationMover();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	virtual std::string get_name() const;

	/// @brief Apply Rotation
	virtual void apply( core::pose::Pose & pose );

private:

	// Store new normal
	core::Real rot_mag_;

	// Store jump num
	core::SSize rb_jump_;

};

/// @brief Random Position Translation Move
class RandomPositionTranslationMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Construct a Default Random Position Mover
	RandomPositionTranslationMover();

	/// @brief Custom Constructor
	/// @details Specify a magnitude of translation
	RandomPositionTranslationMover(
		core::Real trans_mag,
		core::SSize rb_jump = 2
	);

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	RandomPositionTranslationMover( RandomPositionTranslationMover const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	RandomPositionTranslationMover &
	operator=( RandomPositionTranslationMover const & src );

	/// @brief Destructor
	~RandomPositionTranslationMover();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	virtual std::string get_name() const;

	/// @brief Apply Translation to membrane position
	/// @brief Translate membrane position to new center
	virtual void apply( core::pose::Pose & pose );

private:

	// Store trans mag
	core::Real trans_mag_;

	// Store jump num
	core::SSize rb_jump_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_RandomMembranePositionMover_hh
