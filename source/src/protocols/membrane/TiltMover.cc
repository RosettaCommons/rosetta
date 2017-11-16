// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/TiltMover.cc
/// @brief      Tilts a protein in the membrane (Rosetta Scripts Hook)
/// @details Tilts a span, protein or part of a pose in the membrane,
///    depending on the jump number. The tilt axis is the axis
///    perpendicular to the axis connecting the embedding centers of the
///    two partners;
///    BEWARE: CANNOT USE MEMBRANE JUMP AS JUMP NUMBER!!!
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TiltMover_cc
#define INCLUDED_protocols_membrane_TiltMover_cc

// Unit Headers
#include <protocols/membrane/TiltMover.hh>
#include <protocols/membrane/TiltMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/util.hh>

// Package Headers
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
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.TiltMover" );

namespace protocols {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Defaults: jump = 1, angle = random, axis =
/// axis perpendicular to axis connecting protein embedding centers
TiltMover::TiltMover() :
	protocols::moves::Mover()
{
	set_defaults();
	register_options();
}

/// @brief Custom Constructor
/// @details User can specify jump number
TiltMover::TiltMover( core::Size jump_num ) :
	protocols::moves::Mover()
{
	set_defaults();
	register_options();

	jump_num_ = jump_num;
}

/// @brief Custom constructor
/// @details User can specify jump number and angle
TiltMover::TiltMover( core::Size jump_num, core::Real angle ) :
	protocols::moves::Mover()
{
	set_defaults();
	register_options();

	jump_num_ = jump_num;
	angle_ = angle;
	random_angle_ = false;
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
TiltMover::TiltMover( TiltMover const & ) = default;

/// @brief Assignment Operator
TiltMover & TiltMover::operator = ( TiltMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new TiltMover( *this ) );
}

/// @brief Destructor
TiltMover::~TiltMover() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
TiltMover::clone() const {
	return ( protocols::moves::MoverOP( new TiltMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
TiltMover::fresh_instance() const {
	return protocols::moves::MoverOP( new TiltMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
TiltMover::parse_my_tag(
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
TiltMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TiltMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
TiltMoverCreator::keyname() const {
	return TiltMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
TiltMoverCreator::mover_name() {
	return "TiltMover";
}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (TiltMover)
std::string
TiltMover::get_name() const {
	return "TiltMover";
}

/// @brief Flip the downstream partner in the membrane
void TiltMover::apply( core::pose::Pose & pose ) {

	using namespace core;
	using namespace numeric;
	using namespace core::conformation::membrane;
	using namespace protocols::rigid;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;

	TR << "Tilting along a jump in the membrane..." << std::endl;

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// if random angle, set it to random
	if ( random_angle_ == true ) {
		angle_ = numeric::random::random_range( -45, 45 );
		TR << "angle is random" << std::endl;
	}

	TR << "angle is " << angle_ << std::endl;

	// compute downstream empedding
	SpanningTopologyOP topo = pose.conformation().membrane_info()->spanning_topology();
	SpanningTopologyOP topo_up_( new SpanningTopology() );
	SpanningTopologyOP topo_down_( new SpanningTopology() );

	// split_topology_by_jump_noshift
	split_topology_by_jump_noshift( pose, jump_num_, topo, topo_up_, topo_down_ );

	// compute embedding for partners (compute structure-based embedding with split topologies)
	EmbeddingDefOP emb_up( compute_structure_based_embedding( pose, *topo_up_ ) );
	EmbeddingDefOP emb_down( compute_structure_based_embedding( pose, *topo_down_ ) );

	// compute tilt axis
	core::Vector emb_cnt_vector = emb_down->center() - emb_up->center();
	core::Vector mem_normal = pose.conformation().membrane_info()->membrane_normal(pose.conformation());
	core::Vector tangent_axis = cross( emb_cnt_vector, mem_normal );

	// tilt the downstream partner
	RigidBodyDeterministicSpinMoverOP tilt( new RigidBodyDeterministicSpinMover(
		jump_num_, tangent_axis, emb_down->center(), angle_ ) );
	tilt->apply( pose );

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

/// @brief Set Random tilt angle between -20 and 20 degrees
void TiltMover::set_random_membrane_tilt_angle() {
	random_angle_ = true;
} // set random membrane flip angle

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void TiltMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::mp::setup::spanfiles );

}

/// @brief Set default values
void TiltMover::set_defaults() {

	jump_num_ = 1;
	angle_ = 10;
	random_angle_ = true;

}// set_defaults


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TiltMover_cc
