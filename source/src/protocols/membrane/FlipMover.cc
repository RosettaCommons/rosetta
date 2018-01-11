// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/FlipMoverCreator.hh
/// @brief      Flips a span or protein in the membrane (Rosetta Scripts Hook)
/// @details Flips a span, protein or part of a pose in the membrane,
///    depending on the jump number.
///    ONLY FOR FIXED MEMBRANE AND FLEXIBLE PROTEIN
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
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/util.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.membrane.FlipMover" );

namespace protocols {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Defaults: jump = membrane jump, angle = 180 deg,
///   axis = x-axis
FlipMover::FlipMover() : protocols::moves::Mover()
{
	set_defaults();
	register_options();
}

/// @brief Custom Constructor
/// @details User can specify jump number
FlipMover::FlipMover( core::Size jump_num )
{
	set_defaults();
	register_options();

	jump_num_ = jump_num;
}

/// @brief Custom constructor
/// @details User can specify jump number and rotation axis
FlipMover::FlipMover( core::Size jump_num, core::Vector axis )
{
	set_defaults();
	register_options();

	jump_num_ = jump_num;
	axis_ = axis;
}

/// @brief Custom constructor
/// @details User can specify jump number and angle
FlipMover::FlipMover( core::Size jump_num, core::Real angle )
{
	set_defaults();
	register_options();

	jump_num_ = jump_num;
	angle_ = angle;
}

/// @brief Custom constructor
/// @details User can specify jump number and rotation axis
FlipMover::FlipMover( core::Size jump_num, core::Vector axis, core::Real angle )
{
	set_defaults();
	register_options();

	jump_num_ = jump_num;
	axis_ = axis;
	angle_ = angle;
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
FlipMover::FlipMover( FlipMover const & ) = default;

/// @brief Assignment Operator
FlipMover & FlipMover::operator = ( FlipMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new FlipMover( *this ) );
}

/// @brief Destructor
FlipMover::~FlipMover() = default;

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
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FlipMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FlipMover() );
// XRW TEMP }

/// @brief Return the Name of this mover (as seen by Rscripts)
// XRW TEMP std::string
// XRW TEMP FlipMoverCreator::keyname() const {
// XRW TEMP  return FlipMover::mover_name();
// XRW TEMP }

/// @brief Mover name for Rosetta Scripts
// XRW TEMP std::string
// XRW TEMP FlipMover::mover_name() {
// XRW TEMP  return "FlipMover";
// XRW TEMP }


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (FlipMover)
// XRW TEMP std::string
// XRW TEMP FlipMover::get_name() const {
// XRW TEMP  return "FlipMover";
// XRW TEMP }

/// @brief Flip the downstream partner in the membrane
void FlipMover::apply( core::pose::Pose & pose ) {

	using namespace numeric;
	using namespace core::pose;
	using namespace core::conformation::membrane;
	using namespace protocols::rigid;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;

	TR << "Flipping along a jump in the membrane..." << std::endl;

	// if random angle, set it to random
	if ( random_angle_ == true ) {
		if ( max_angle_dev_ == 0 ) {
			angle_ = numeric::random::random_range( 135, 225 );
		} else {
			angle_ = numeric::random::random_range( 180 - max_angle_dev_, 180 + max_angle_dev_ );
		}
	}

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// reorder foldtree
	Size mem_rsd = pose.membrane_info()->membrane_rsd_num();
	core::kinematics::FoldTree foldtree = pose.fold_tree();
	foldtree.reorder( mem_rsd );
	pose.fold_tree( foldtree );
	TR << "foldtree reordered" << std::endl;
	pose.fold_tree().show(TR);

	TR << "upstream jump res: " << foldtree.upstream_jump_residue( jump_num_ ) << std::endl;
	TR << "downstream jump res: " << foldtree.downstream_jump_residue( jump_num_ ) << std::endl;
	TR << "membrane res: " << pose.membrane_info()->membrane_rsd_num() << std::endl;

	// if jump is not membrane jump
	if ( foldtree.upstream_jump_residue( jump_num_ ) != mem_rsd &&
			foldtree.downstream_jump_residue( jump_num_ ) != mem_rsd ) {

		TR << "jump is not membrane jump" << std::endl;

		// split pose and topology by jump
		Pose pose_up, pose_down;
		SpanningTopologyOP topo = pose.membrane_info()->spanning_topology();
		SpanningTopology topo_up, topo_down;

		// call function
		split_topology_by_jump( pose, jump_num_, *topo, pose_up, pose_down, topo_up, topo_down );

		// get embedding centers of both partners
		EmbeddingDefOP emb_pose_up = compute_structure_based_embedding( pose_up, topo_up );
		EmbeddingDefOP emb_pose_down = compute_structure_based_embedding( pose_down, topo_down );

		// get rotation axis, this is the axis between the embedding centers from both partners
		if ( axis_.length() == 0 ) {
			axis_ = emb_pose_up->center() - emb_pose_down->center();
		}

		// rotation center is embedding center of downstream pose
		core::Vector rot_center = emb_pose_down->center();

		TR << "jump_num_: " << jump_num_ << std::endl;
		TR << "axis_: " << axis_.to_string() << std::endl;
		TR << "rot_center: " << rot_center.to_string() << std::endl;
		TR << "angle: " << angle_ << std::endl;

		// do rotation
		RigidBodyDeterministicSpinMoverOP rb_flip ( new RigidBodyDeterministicSpinMover( jump_num_, axis_, rot_center, angle_ ) );
		rb_flip->apply( pose );
	} else {
		// if jump is somehow connected to the membrane

		TR << "jump is either zero or membrane jump" << std::endl;

		jump_num_ = pose.membrane_info()->membrane_jump();

		// axis is x-axis
		if ( axis_.length() == 0 ) {
			axis_.assign(1, 0, 0);
		}

		// axis is axis between membrane center and pose embedding center
		EmbeddingDefOP embedding( compute_structure_based_embedding( pose ) );
		core::Vector rot_center( embedding->center() );

		TR << "jump_num_: " << jump_num_ << std::endl;
		TR << "axis_: " << axis_.to_string() << std::endl;
		TR << "rot_center: " << rot_center.to_string() << std::endl;
		TR << "angle: " << angle_ << std::endl;

		// do rotation, rot center is pose embedding center
		RigidBodyDeterministicSpinMoverOP rb_flip ( new RigidBodyDeterministicSpinMover( jump_num_, axis_, rot_center, angle_ ) );
		rb_flip->apply( pose );
	}

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

/// @brief Set Random flip angle between 160 and 200 degrees to keep
///   protein oriented in the membrane correctly
void FlipMover::set_random_membrane_flip_angle() {
	random_angle_ = true;
} // set random membrane flip angle

/// @brief Set angle range
/// @details Maximum angle deviation from 180 degrees
void FlipMover::set_range( core::Real max_angle_dev ) {
	max_angle_dev_ = max_angle_dev;
}

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void FlipMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::mp::setup::spanfiles );

}

/// @brief Set default values
void FlipMover::set_defaults() {

	jump_num_ = 1;

	// default rotation axis is x-axis
	axis_.assign( 0, 0, 0 );
	angle_ = 180;

	random_angle_ = false;

	max_angle_dev_ = 0;

}// set_defaults

std::string FlipMover::get_name() const {
	return mover_name();
}

std::string FlipMover::mover_name() {
	return "FlipMover";
}

void FlipMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Flips pose in the membrane", attlist );
}

std::string FlipMoverCreator::keyname() const {
	return FlipMover::mover_name();
}

protocols::moves::MoverOP
FlipMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new FlipMover );
}

void FlipMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FlipMover::provide_xml_schema( xsd );
}



} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_FlipMover_cc
