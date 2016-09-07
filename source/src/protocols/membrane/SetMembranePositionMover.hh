// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/SetMembranePositionMover.hh
///
/// @brief  Sets the membrane position normal and center
/// @details Sets the membrane position normal and center
///    CAUTION: ONLY FOR FLEXIBLE MEMBRANE AND FIXED PROTEIN!!!
///    Last Modified: 6/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_SetMembranePositionMover_hh
#define INCLUDED_protocols_membrane_SetMembranePositionMover_hh

// Unit Headers
#include <protocols/membrane/SetMembranePositionMover.fwd.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/UniformPositionMover.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>

namespace protocols {
namespace membrane {

/// @brief Membrane Position Translation-Rotation Mover
class SetMembranePositionMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Construct a Default Membrane Position Mover
	SetMembranePositionMover();

	/// @brief Custom Constructor
	/// @details Specify a new membrane center and normal
	/// to move this position to
	SetMembranePositionMover( core::Vector center, core::Vector normal );

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	SetMembranePositionMover( SetMembranePositionMover const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	SetMembranePositionMover &
	operator=( SetMembranePositionMover const & src );

	/// @brief Destructor
	~SetMembranePositionMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	std::string get_name() const override;

	/// @brief Apply Rotation/Translation to Membrane
	/// @brief Translate the membrane position in this pose
	/// to the new center position, and rotate to new normal
	void apply( core::pose::Pose & pose ) override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

private:

	// Store new normal/center pair
	core::Vector center_;
	core::Vector normal_;

};

/// @brief Membrane Position Rotation Move
/// @details Rotate the orientation of the membrane position to a new
/// normal position. Applies rotation to jump
class SetMembraneNormalMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Construct a Default Membrane Position Mover
	SetMembraneNormalMover();

	/// @brief Custom Constructor
	/// @details Specify a new normal to rotate membranes to
	/// to move this position to
	SetMembraneNormalMover( core::Vector normal );

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	SetMembraneNormalMover( SetMembraneNormalMover const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	SetMembraneNormalMover &
	operator=( SetMembraneNormalMover const & src );

	/// @brief Destructor
	~SetMembraneNormalMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	std::string get_name() const override;

	/// @brief Apply Rotation to Membrane
	/// @brief Rotate the membrane to the new normal position
	void apply( core::pose::Pose & pose ) override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

private:

	// Store new normal
	core::Vector normal_;

};

/// @brief Membrane Position Translation Move
/// @details Translate the center of the membrane stub ot the specified position
class SetMembraneCenterMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Construct a Default Membrane Position Mover
	SetMembraneCenterMover();

	/// @brief Custom Constructor
	/// @details Specify a new center position to translate this stub to
	SetMembraneCenterMover( core::Vector center );

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	SetMembraneCenterMover( SetMembraneCenterMover const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	SetMembraneCenterMover &
	operator=( SetMembraneCenterMover const & src );

	/// @brief Destructor
	~SetMembraneCenterMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	std::string get_name() const override;

	/// @brief Apply Translation to membrane position
	/// @brief Translate membrane position to new center
	void apply( core::pose::Pose & pose ) override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

private:

	// Store new center
	core::Vector center_;

};

} // membrane
} // protocols

#endif // #ifndef INCLUDED_protocols_membrane_SetMembranePositionMover_hh
