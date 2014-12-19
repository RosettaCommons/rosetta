// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/FlipMoverCreator.hh
/// @brief      Flips a span or protein in the membrane (Rosetta Scripts Hook)
/// @details	Flips a span, protein or part of a pose in the membrane,
///				depending on the jump number.
///				ONLY FOR FIXED MEMBRANE AND FLEXIBLE PROTEIN
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_FlipMover_cc
#define INCLUDED_protocols_membrane_FlipMover_cc

// Unit Headers
#include <protocols/membrane/FlipMover.hh> 
#include <protocols/membrane/FlipMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/util.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh> 
#include <core/pose/util.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <core/conformation/membrane/types.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.FlipMover" );

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;
		
/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Defaults: jump = membrane jump, angle = 180 deg, axis =
///			 axis between COMs projected into the membrane plane
FlipMover::FlipMover() : protocols::moves::Mover()
{
	set_defaults();
	register_options();
}

///// @brief Custom Constructor
///// @details User can specify jump number
//FlipMover::FlipMover( Size jump_num )
//{
//	set_defaults();
//	register_options();
//	
//	jump_num_ = jump_num;
//}
//
///// @brief Custom constructor
///// @details User can specify jump number and rotation axis
//FlipMover::FlipMover( Size jump_num, Vector axis )
//{
//	set_defaults();
//	register_options();
//	
//	jump_num_ = jump_num;
//	axis_ = axis;
//}
//
///// @brief Custom constructor
///// @details User can specify jump number and angle
//FlipMover::FlipMover( Size jump_num, Real angle )
//{
//	set_defaults();
//	register_options();
//	
//	jump_num_ = jump_num;
//	angle_ = angle;
//}
//
///// @brief Custom constructor
///// @details User can specify jump number and rotation axis
//FlipMover::FlipMover( Size jump_num, Vector axis, Real angle )
//{
//	set_defaults();
//	register_options();
//	
//	jump_num_ = jump_num;
//	axis_ = axis;
//	angle_ = angle;
//}
//	
///// @brief Copy Constructor
///// @details Create a deep copy of this mover
//FlipMover::FlipMover( FlipMover const & src ) : protocols::moves::Mover( src ),
//	jump_num_( src.jump_num_ ),
//	axis_( src.axis_ ),
//	angle_( src.angle_ )
//{}

/// @brief Assignment Operator
FlipMover & FlipMover::operator = ( FlipMover const & src ) {
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
		
	// Otherwise, create a new object
	return *( new FlipMover( *this ) );
}

/// @brief Destructor
FlipMover::~FlipMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
FlipMover::clone() const {
	return ( protocols::moves::MoverOP( new FlipMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
FlipMover::fresh_instance() const {
	return protocols::moves::MoverOP( new FlipMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
FlipMover::parse_my_tag(
	 utility::tag::TagCOP /*tag*/,
	 basic::datacache::DataMap &,
	 protocols::filters::Filters_map const &,
	 protocols::moves::Movers_map const &,
	 core::pose::Pose const &
	 ) {

	// TODO: implement this
	
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
FlipMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new FlipMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
FlipMoverCreator::keyname() const {
	return FlipMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
FlipMoverCreator::mover_name() {
	return "FlipMover";
}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (FlipMover)
std::string
FlipMover::get_name() const {
	return "FlipMover";
}

/// @brief Flip the downstream partner in the membrane
void FlipMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::conformation::membrane;
	using namespace protocols::rigid;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;
	
	TR << "Flipping along a jump in the membrane..." << std::endl;

	// initialize jump
	if ( jump_num_ == 0 ) {
		jump_num_ = pose.conformation().membrane_info()->membrane_jump();
	}

	// reorder foldtree
	core::kinematics::FoldTree foldtree = pose.fold_tree();
	foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
	pose.fold_tree( foldtree );
	TR << "foldtree reordered" << std::endl;
	pose.fold_tree().show(std::cout);

	// set axis to x axis
	axis_.assign( 1, 0, 0 );
	
	// get the membrane axis
//	if ( axis_.length() == 0 ) {
//		axis_ = membrane_axis( pose, jump_num_ );
//	}
//
//	// split pose into subposes by jump
//	Pose upstream_pose, downsteam_pose;
//	partition_pose_by_jump( pose, jump_num_, upstream_pose, downsteam_pose );
//	TR << "up pose: " << upstream_pose.total_residue() << std::endl;
//	TR << "down pose: " << downsteam_pose.total_residue() << std::endl;

	// rotation center is total embedding center
	EmbeddingDefOP embedding( compute_structure_based_membrane_position( pose ) );
	Vector rot_center = embedding->center();
	
	// TODO:
	// use partition_pose_by_jump in core/pose/util
	// compute_structure_based_embedding from subpose
	// rotation center is center of total embedding

	

//	void
//	partition_pose_by_jump(
//						   pose::Pose const & src,
//						   int const jump_number,
//						   pose::Pose & partner1,
//						   pose::Pose & partner2
//						   );
	
	
//	// get rotation center
//	core::kinematics::Stub downstream_stub = pose.conformation().downstream_jump_stub( jump_num_ );
//	Vector rot_center = downstream_stub.center();
//	rot_center.z( 0 );

//	Vector dummy_up, dummy_down;
//	protocols::geometry::centroids_by_jump(pose, rb_jump_, dummy_up, dummy_down);
//	rot_center = dummy_down;

	// do rotation
	RigidBodyDeterministicSpinMoverOP rb_flip ( new RigidBodyDeterministicSpinMover( jump_num_, axis_, rot_center, angle_ ) );
	rb_flip->apply( pose );
	
}// apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// membrane_new, seutp options: center, normal, spanfile and
void FlipMover::register_options() {
	
	using namespace basic::options;
	option.add_relevant( OptionKeys::membrane_new::setup::spanfiles );
	
}

/// @brief Set default values
void FlipMover::set_defaults() {
	
	jump_num_ = 0;
	axis_.assign( 0, 0, 0 );
	angle_ = 180;
	
}// set_defaults


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_FlipMover_cc
