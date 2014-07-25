// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/MembranePositionRTMover.cc
///
/// @brief		Membrane Position Rotation/Translation Mover
///	@details	Apply a uniform rigid translation & rotation of the
///				membrane defined by the center/normal coordinates
///				stored in the membrane.
///				Last Modified: 6/28/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <protocols/membrane/MembranePositionRTMover.hh>
#include <protocols/membrane/MembranePositionRTMoverCreator.hh>

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

static basic::Tracer TR( "protocols.membrane.MembranePositionRTMOver" );

namespace protocols {
namespace membrane {

using namespace core;
using namespace protocols::membrane;
using namespace protocols::simple_moves;

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
MembranePositionRTMover::MembranePositionRTMover() :
	Mover(),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane center and normal
///	to move this position to
MembranePositionRTMover::MembranePositionRTMover( Vector center, Vector normal ) :
	 Mover(),
	 center_( center ),
	 normal_( normal )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
MembranePositionRTMover::MembranePositionRTMover( MembranePositionRTMover const & src ) :
	Mover( src ),
	center_( src.center_ ),
	normal_( src.normal_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
MembranePositionRTMover &
MembranePositionRTMover::operator=( MembranePositionRTMover const & src )
{
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new MembranePositionRTMover( *this ) );
	
}

/// @brief Destructor
MembranePositionRTMover::~MembranePositionRTMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
MembranePositionRTMover::get_name() const {
	return "MembranePositionRTMover";
}

/// @brief Apply Rotation/Translation to Membrane
/// @brief Translate the membrane position in this pose
/// to the new center position, and rotate to new normal
void
MembranePositionRTMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	
	// Check the pose is ==a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}
	
	// Check the membrane fold tree is reasonable
	if (! pose.conformation().membrane()->check_membrane_fold_tree( pose.fold_tree() ) ) {
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
MembranePositionRTMover::clone() const {
	return ( new MembranePositionRTMover( *this ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MembranePositionRTMover::fresh_instance() const {
	return new MembranePositionRTMover();
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MembranePositionRTMover::parse_my_tag(
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
MembranePositionRTMoverCreator::create_mover() const {
	return new MembranePositionRTMover;
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MembranePositionRTMoverCreator::keyname() const {
	return MembranePositionRTMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MembranePositionRTMoverCreator::mover_name() {
	return "MembranePositionRTMoverr";
}

// Membrane Position Rotation Mover /////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
MembranePositionRotationMover::MembranePositionRotationMover() :
	Mover(),
	normal_( 0, 0, 1 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane normal
///	to move this position to
MembranePositionRotationMover::MembranePositionRotationMover( Vector normal ) :
	Mover(),
	normal_( normal )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
MembranePositionRotationMover::MembranePositionRotationMover( MembranePositionRotationMover const & src ) :
	Mover( src ),
	normal_( src.normal_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
MembranePositionRotationMover &
MembranePositionRotationMover::operator=( MembranePositionRotationMover const & src )
{
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new MembranePositionRotationMover( *this ) );
	
}

/// @brief Destructor
MembranePositionRotationMover::~MembranePositionRotationMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
MembranePositionRotationMover::get_name() const {
	return "MembranePositionRotationMover";
}

/// @brief Apply Rotation move to Membrane
/// @brief Rotate membrane position to new normal
void
MembranePositionRotationMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	
	// Check the pose is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}
	
	// Check the membrane fold tree is reasonable
	if (! pose.conformation().membrane()->check_membrane_fold_tree( pose.fold_tree() ) ) {
		utility_exit_with_message( "Cannot apply membrane move with unreasonable membrane fold tree" );
	}
	
	// Compute Rotation Axis - CrossProd between Current & New Normal axis
	Vector current_normal( pose.conformation().membrane_normal() );
	Vector current_center( pose.conformation().membrane_center() );
	
	pose.conformation().update_membrane_position( current_center, current_center + normal_ );

}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MembranePositionRotationMover::clone() const {
	return ( new MembranePositionRotationMover( *this ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MembranePositionRotationMover::fresh_instance() const {
	return new MembranePositionRotationMover();
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MembranePositionRotationMover::parse_my_tag(
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
MembranePositionRotationMoverCreator::create_mover() const {
	return new MembranePositionRotationMover;
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MembranePositionRotationMoverCreator::keyname() const {
	return MembranePositionRotationMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MembranePositionRotationMoverCreator::mover_name() {
	return "MembranePositionRotationMoverr";
}

// Membrane Position Rotation Mover /////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
MembranePositionTranslationMover::MembranePositionTranslationMover() :
	Mover(),
	center_( 0, 0, 0 )
{}

/// @brief Custom Constructor
/// @details Specify a new membrane center to move to
MembranePositionTranslationMover::MembranePositionTranslationMover( Vector center ) :
	Mover(),
	center_( center )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
MembranePositionTranslationMover::MembranePositionTranslationMover( MembranePositionTranslationMover const & src ) :
	Mover( src ),
	center_( src.center_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
MembranePositionTranslationMover &
MembranePositionTranslationMover::operator=( MembranePositionTranslationMover const & src )
{
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new MembranePositionTranslationMover( *this ) );
	
}

/// @brief Destructor
MembranePositionTranslationMover::~MembranePositionTranslationMover() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
MembranePositionTranslationMover::get_name() const {
	return "MembranePositionTranslationMover";
}

/// @brief Apply Rotation/Translation to Membrane
/// @brief Translate the membrane position in this pose
/// to the new center position, and rotate to new normal
void
MembranePositionTranslationMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	
	// Check the pose is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}
	
	// Check the membrane fold tree is reasonable
	if (! pose.conformation().membrane()->check_membrane_fold_tree( pose.fold_tree() ) ) {
		utility_exit_with_message( "Cannot apply membrane move with unreasonable membrane fold tree" );
	}
	
	// Get current normal
	Vector current_normal( pose.conformation().membrane_normal() );
	
	// Apply translation
	pose.conformation().update_membrane_position( center_, current_normal + center_ );

}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MembranePositionTranslationMover::clone() const {
	return ( new MembranePositionTranslationMover( *this ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MembranePositionTranslationMover::fresh_instance() const {
	return new MembranePositionTranslationMover();
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MembranePositionTranslationMover::parse_my_tag(
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
MembranePositionTranslationMoverCreator::create_mover() const {
	return new MembranePositionTranslationMover;
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MembranePositionTranslationMoverCreator::keyname() const {
	return MembranePositionTranslationMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MembranePositionTranslationMoverCreator::mover_name() {
	return "MembranePositionTranslationMoverr";
}

} // membrane
} // protocols

