// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/membrane/SetMembranePositionMover.cc
///
/// @brief  Sets the membrane position normal and center
///    Last Modified: 7/4/15
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/SetMembranePositionMoverCreator.hh>

// Project Headers
#include <protocols/membrane/util.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/simple_moves/UniformPositionMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Package Headers
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <utility/tag/Tag.hh>

static thread_local basic::Tracer TR( "protocols.membrane.SetMembranePositionMover" );

namespace protocols {
namespace membrane {

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
SetMembranePositionMover::SetMembranePositionMover() :
	protocols::moves::Mover(),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane center and normal
/// to move this position to
SetMembranePositionMover::SetMembranePositionMover( core::Vector center, core::Vector normal ) :
	protocols::moves::Mover(),
	center_( center ),
	normal_( normal )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
SetMembranePositionMover::SetMembranePositionMover( SetMembranePositionMover const & src ) :
	protocols::moves::Mover( src ),
	center_( src.center_ ),
	normal_( src.normal_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
SetMembranePositionMover &
SetMembranePositionMover::operator=( SetMembranePositionMover const & src )
{

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new SetMembranePositionMover( *this ) );

}

/// @brief Destructor
SetMembranePositionMover::~SetMembranePositionMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
SetMembranePositionMover::get_name() const {
	return "SetMembranePositionMover";
}

/// @brief Apply Rotation/Translation to Membrane
/// @brief Translate the membrane position in this pose
/// to the new center position, and rotate to new normal
void
SetMembranePositionMover::apply( core::pose::Pose & pose ) {

	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;

	TR << "Calling SetMembranePositionMover" << std::endl;

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// Check the pose is ==a membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}

	// Update coordinates
	pose.conformation().update_membrane_position( center_, normal_ );

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

} // apply

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
SetMembranePositionMover::clone() const {
	return ( protocols::moves::MoverOP( new SetMembranePositionMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
SetMembranePositionMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SetMembranePositionMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
SetMembranePositionMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	read_center_normal_from_tag( center_, normal_, tag );
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
SetMembranePositionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetMembranePositionMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
SetMembranePositionMoverCreator::keyname() const {
	return SetMembranePositionMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
SetMembranePositionMoverCreator::mover_name() {
	return "SetMembranePositionMover";
}

// SetMembraneNormalMover /////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
SetMembraneNormalMover::SetMembraneNormalMover() :
	protocols::moves::Mover(),
	normal_( 0, 0, 1 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane normal
/// to move this position to
SetMembraneNormalMover::SetMembraneNormalMover( core::Vector normal ) :
	protocols::moves::Mover(),
	normal_( normal )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
SetMembraneNormalMover::SetMembraneNormalMover( SetMembraneNormalMover const & src ) :
	protocols::moves::Mover( src ),
	normal_( src.normal_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
SetMembraneNormalMover &
SetMembraneNormalMover::operator=( SetMembraneNormalMover const & src )
{

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new SetMembraneNormalMover( *this ) );

}

/// @brief Destructor
SetMembraneNormalMover::~SetMembraneNormalMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
SetMembraneNormalMover::get_name() const {
	return "SetMembraneNormalMover";
}

/// @brief Apply Rotation move to Membrane
/// @brief Rotate membrane position to new normal
void
SetMembraneNormalMover::apply( core::pose::Pose & pose ) {

	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;

	// Check the pose is a membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}

	// Compute Rotation Axis - CrossProd between Current & New Normal axis
	core::Vector current_normal( pose.conformation().membrane_info()->membrane_normal() );
	core::Vector current_center( pose.conformation().membrane_info()->membrane_center() );

	pose.conformation().update_membrane_position( current_center, normal_ );

}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
SetMembraneNormalMover::clone() const {
	return ( protocols::moves::MoverOP( new SetMembraneNormalMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
SetMembraneNormalMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SetMembraneNormalMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
SetMembraneNormalMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	core::Vector center_dummy(0,0,0);
	read_center_normal_from_tag( center_dummy, normal_, tag );

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
SetMembraneNormalMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetMembraneNormalMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
SetMembraneNormalMoverCreator::keyname() const {
	return SetMembraneNormalMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
SetMembraneNormalMoverCreator::mover_name() {
	return "SetMembraneNormalMover";
}

// SetMembraneCenterMover /////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
SetMembraneCenterMover::SetMembraneCenterMover() :
	protocols::moves::Mover(),
	center_( 0, 0, 0 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane center to move to
SetMembraneCenterMover::SetMembraneCenterMover( core::Vector center ) :
	protocols::moves::Mover(),
	center_( center )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
SetMembraneCenterMover::SetMembraneCenterMover( SetMembraneCenterMover const & src ) :
	protocols::moves::Mover( src ),
	center_( src.center_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
SetMembraneCenterMover &
SetMembraneCenterMover::operator=( SetMembraneCenterMover const & src )
{

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new SetMembraneCenterMover( *this ) );

}

/// @brief Destructor
SetMembraneCenterMover::~SetMembraneCenterMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
SetMembraneCenterMover::get_name() const {
	return "SetMembraneCenterMover";
}

/// @brief Apply Rotation/Translation to Membrane
/// @brief Translate the membrane position in this pose
/// to the new center position, and rotate to new normal
void
SetMembraneCenterMover::apply( core::pose::Pose & pose ) {

	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;

	// Check the pose is a membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}

	// Get current normal
	core::Vector current_normal( pose.conformation().membrane_info()->membrane_normal() );

	// Apply translation
	pose.conformation().update_membrane_position( center_, current_normal );

}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
SetMembraneCenterMover::clone() const {
	return ( protocols::moves::MoverOP( new SetMembraneCenterMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
SetMembraneCenterMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SetMembraneCenterMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
SetMembraneCenterMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	core::Vector normal_dummy(0,0,0);
	read_center_normal_from_tag( center_, normal_dummy, tag );

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
SetMembraneCenterMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetMembraneCenterMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
SetMembraneCenterMoverCreator::keyname() const {
	return SetMembraneCenterMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
SetMembraneCenterMoverCreator::mover_name() {
	return "SetMembraneCenterMoverr";
}

} // membrane
} // protocols

