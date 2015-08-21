// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/SpinAroundPartnerMover.cc
/// @brief      Spins the downstream partner around the upstream partner
/// @details Spins the downstream partner around the upstream partner in the
///    membrane to probe all kinds of interfaces. Both embedding normals
///    are approximately conserved, i.e. the partners aren't flipped
///    in the membrane.
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_SpinAroundPartnerMover_cc
#define INCLUDED_protocols_membrane_SpinAroundPartnerMover_cc

// Unit Headers
#include <protocols/membrane/SpinAroundPartnerMover.hh>
#include <protocols/membrane/SpinAroundPartnerMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
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

static basic::Tracer TR( "protocols.membrane.SpinAroundPartnerMover" );

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
/// @details Defaults: jump = 1, sampling range = 100
///    Sampling range of 100 means that both x and y are sampled from
///    -100 to +100 before calling DockingSlideIntoContact
SpinAroundPartnerMover::SpinAroundPartnerMover() : protocols::moves::Mover()
{
	set_defaults();
	register_options();
}

/// @brief Custom Constructor
/// @details User can specify jump number
SpinAroundPartnerMover::SpinAroundPartnerMover( Size jump_num )
{
	set_defaults();
	register_options();

	jump_ = jump_num;
}

/// @brief Custom constructor
/// @details User can specify jump number and sampling range
SpinAroundPartnerMover::SpinAroundPartnerMover( Size jump_num, Size range )
{
	set_defaults();
	register_options();

	jump_ = jump_num;
	rand_range_ = true;
	range_ = range;
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
SpinAroundPartnerMover::SpinAroundPartnerMover( SpinAroundPartnerMover const & src ) : protocols::moves::Mover( src ),
	jump_( src.jump_ ),
	rand_range_( src.rand_range_ ),
	range_( src.range_ ),
	x_( src.x_ ),
	y_( src.y_ )
{}

/// @brief Assignment Operator
SpinAroundPartnerMover & SpinAroundPartnerMover::operator = ( SpinAroundPartnerMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new SpinAroundPartnerMover( *this ) );
}

/// @brief Destructor
SpinAroundPartnerMover::~SpinAroundPartnerMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
SpinAroundPartnerMover::clone() const {
	return ( protocols::moves::MoverOP( new SpinAroundPartnerMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
SpinAroundPartnerMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SpinAroundPartnerMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
SpinAroundPartnerMover::parse_my_tag(
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
SpinAroundPartnerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SpinAroundPartnerMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
SpinAroundPartnerMoverCreator::keyname() const {
	return SpinAroundPartnerMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
SpinAroundPartnerMoverCreator::mover_name() {
	return "SpinAroundPartnerMover";
}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (SpinAroundPartnerMover)
std::string
SpinAroundPartnerMover::get_name() const {
	return "SpinAroundPartnerMover";
}

/// @brief Set random range
void SpinAroundPartnerMover::random_range( bool yesno ) {
	rand_range_ = yesno;
}

/// @brief Set x position
void SpinAroundPartnerMover::set_x( Real x ) {
	x_ = x;
	rand_range_ = false;
}

/// @brief Set y position
void SpinAroundPartnerMover::set_y( Real y ) {
	y_ = y;
	rand_range_ = false;
}


/// @brief Flip the downstream partner in the membrane
void SpinAroundPartnerMover::apply( Pose & pose ) {

	using namespace numeric;
	using namespace core::conformation::membrane;
	using namespace protocols::rigid;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;

	TR << "Spinning downstream partner around upstream partner in the membrane..." << std::endl;

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// compute downstream empedding
	SpanningTopologyOP topo = pose.conformation().membrane_info()->spanning_topology();
	SpanningTopologyOP topo_up_( new SpanningTopology() );
	SpanningTopologyOP topo_down_( new SpanningTopology() );

	// split_topology_by_jump_noshift
	split_topology_by_jump_noshift( pose, jump_, topo, topo_up_, topo_down_ );

	// compute embedding for partners (compute structure-based embedding with split topologies)
	EmbeddingDefOP emb_up( compute_structure_based_embedding( pose, *topo_up_ ) );
	EmbeddingDefOP emb_down( compute_structure_based_embedding( pose, *topo_down_ ) );

	// get random x and y positions
	if ( rand_range_ == true ) {

		x_ = static_cast< core::Real > ( numeric::random::random_range( 0, range_ ) );
		y_ = static_cast< core::Real > ( numeric::random::random_range( 0, range_ ) );
		x_ -= ( range_ / 2 );
		y_ -= ( range_ / 2 );
	}

	// set new embedding vector of downstream partner
	core::Vector new_emb_center( x_, y_, 0 );
	TR << "translation: new emb center: " << new_emb_center.to_string() << std::endl;
	core::Vector trans_vector = new_emb_center - emb_down->center();
	TR << "translation: vector: " << trans_vector.to_string() << std::endl;
	TranslationMoverOP trans( new TranslationMover( trans_vector, jump_ ) );
	trans->apply( pose );

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void SpinAroundPartnerMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::mp::setup::spanfiles );

}

/// @brief Set default values
void SpinAroundPartnerMover::set_defaults() {

	jump_ = 1;
	rand_range_ = true;
	range_ = 100;
	x_ = 50.0;
	y_ = 50.0;

}// set_defaults


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_SpinAroundPartnerMover_cc
