// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/SetMembranePositionMover.cc
///
/// @brief		Membrane Position Rotation/Translation Mover
///	@details	Apply a uniform rigid translation & rotation of the
///				membrane defined by the center/normal coordinates
///				stored in the membrane.
///				Last Modified: 6/28/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/SetMembranePositionMoverCreator.hh>

// Project Headers
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

using namespace core;
using namespace protocols::membrane;
using namespace protocols::simple_moves;

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
SetMembranePositionMover::SetMembranePositionMover() :
	Mover(),
	center_( 0.0, 0.0, 0.0 ),
	normal_( 0.0, 0.0, 1.0 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane center and normal
///	to move this position to
SetMembranePositionMover::SetMembranePositionMover( Vector center, Vector normal ) :
	 Mover(),
	 center_( center ),
	 normal_( normal )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
SetMembranePositionMover::SetMembranePositionMover( SetMembranePositionMover const & src ) :
	Mover( src ),
	center_( src.center_ ),
	normal_( src.normal_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
SetMembranePositionMover &
SetMembranePositionMover::operator=( SetMembranePositionMover const & src )
{
	
	// Abort self-assignment.
	if (this == &src) {
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
SetMembranePositionMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	
	// Check the pose is ==a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}
	
	// Check the membrane fold tree is reasonable
	if (! pose.conformation().membrane_info()->check_membrane_fold_tree( pose.fold_tree() ) ) {
		utility_exit_with_message( "Cannot apply membrane move with unreasonable membrane fold tree" );
	}
	
	// Update coordinates
	pose.conformation().update_membrane_position( center_, normal_ + center_ );

}

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
	utility::tag::TagCOP, //tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
	) {
	
	// TODO: implement this
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
	return "SetMembranePositionMoverr";
}

// Membrane Position Rotation Mover /////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
SetMembraneNomalMover::SetMembraneNomalMover() :
	Mover(),
	normal_( 0.0, 0.0, 1.0 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane normal
///	to move this position to
SetMembraneNomalMover::SetMembraneNomalMover( Vector normal ) :
	Mover(),
	normal_( normal )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
SetMembraneNomalMover::SetMembraneNomalMover( SetMembraneNomalMover const & src ) :
	Mover( src ),
	normal_( src.normal_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
SetMembraneNomalMover &
SetMembraneNomalMover::operator=( SetMembraneNomalMover const & src )
{
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new SetMembraneNomalMover( *this ) );
	
}

/// @brief Destructor
SetMembraneNomalMover::~SetMembraneNomalMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
SetMembraneNomalMover::get_name() const {
	return "SetMembraneNomalMover";
}

/// @brief Apply Rotation move to Membrane
/// @brief Rotate membrane position to new normal
void
SetMembraneNomalMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	
	// Check the pose is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}
	
	// Check the membrane fold tree is reasonable
	if (! pose.conformation().membrane_info()->check_membrane_fold_tree( pose.fold_tree() ) ) {
		utility_exit_with_message( "Cannot apply membrane move with unreasonable membrane fold tree" );
	}
	
	// Compute Rotation Axis - CrossProd between Current & New Normal axis
	Vector current_normal( pose.conformation().membrane_info()->membrane_normal() );
	Vector current_center( pose.conformation().membrane_info()->membrane_center() );
	
	pose.conformation().update_membrane_position( current_center, current_center + normal_ );

}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
SetMembraneNomalMover::clone() const {
	return ( protocols::moves::MoverOP( new SetMembraneNomalMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
SetMembraneNomalMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SetMembraneNomalMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
SetMembraneNomalMover::parse_my_tag(
	utility::tag::TagCOP, //tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
	) {
	
	// TODO: implement this
	
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
SetMembraneNomalMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetMembraneNomalMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
SetMembraneNomalMoverCreator::keyname() const {
	return SetMembraneNomalMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
SetMembraneNomalMoverCreator::mover_name() {
	return "SetMembraneNomalMoverr";
}

// Membrane Position Rotation Mover /////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
SetMembraneCenterMover::SetMembraneCenterMover() :
	Mover(),
	center_( 0.0, 0.0, 0.0 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane center to move to
SetMembraneCenterMover::SetMembraneCenterMover( Vector center ) :
	Mover(),
	center_( center )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
SetMembraneCenterMover::SetMembraneCenterMover( SetMembraneCenterMover const & src ) :
	Mover( src ),
	center_( src.center_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
SetMembraneCenterMover &
SetMembraneCenterMover::operator=( SetMembraneCenterMover const & src )
{
	
	// Abort self-assignment.
	if (this == &src) {
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
SetMembraneCenterMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	
	// Check the pose is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}
	
	// Check the membrane fold tree is reasonable
	if (! pose.conformation().membrane_info()->check_membrane_fold_tree( pose.fold_tree() ) ) {
		utility_exit_with_message( "Cannot apply membrane move with unreasonable membrane fold tree" );
	}
	
	// Get current normal
	Vector current_normal( pose.conformation().membrane_info()->membrane_normal() );
	
	// Apply translation
	pose.conformation().update_membrane_position( center_, current_normal + center_ );

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
   utility::tag::TagCOP, // tag,
   basic::datacache::DataMap &,
   protocols::filters::Filters_map const &,
   protocols::moves::Movers_map const &,
   core::pose::Pose const &
   ) {
	
	// TODO: implement this
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

