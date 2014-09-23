// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/RandomOmegaFlipMover.cc
/// @brief RandomOmegaFlipMover methods implemented
/// @author

// unit headers
#include <protocols/simple_moves/RandomOmegaFlipMover.hh>
#include <protocols/simple_moves/RandomOmegaFlipMoverCreator.hh>

// protocols headers
#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/pose/Pose.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/MoveMap.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// numeric headers
#include <numeric/random/random.hh>

// basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.simple_moves.RandomOmegaFlipMover" );

namespace protocols {
namespace simple_moves {


RandomOmegaFlipMover::RandomOmegaFlipMover() :
	Mover("RandomOmegaFlipMover"),
	move_map_( /* 0 */ )
{}

RandomOmegaFlipMover::RandomOmegaFlipMover( core::kinematics::MoveMapOP move_map ) :
	Mover("RandomOmegaFlipMover"),
	move_map_( move_map )
{}

RandomOmegaFlipMover::RandomOmegaFlipMover( RandomOmegaFlipMover const & other ) :
	Mover("RandomOmegaFlipMover"),
	move_map_( core::kinematics::MoveMapOP( new core::kinematics::MoveMap( *other.move_map_ ) ) )
{}

RandomOmegaFlipMover::~RandomOmegaFlipMover(){}

void
RandomOmegaFlipMover::apply( core::pose::Pose & pose )
{
	using namespace basic;
	using namespace core;
	using namespace conformation;
	using namespace kinematics;

	// setup the torsion id list
	setup_torsion_list( pose );

	// empty move maps have zero torsions
	if ( torsion_id_list_.size() ) {

		// randomly select a free torsion
		Size tor_num( numeric::random::rg().random_range( 1, torsion_id_list_.size() ) );

		// calc new value
		Real old_tor( pose.conformation().torsion( torsion_id_list_[ tor_num ] ) );
		Real new_tor( periodic_range( old_tor - 180.0, 360.0 ) );

		//TR << "DEBUG: Setting torsion " << torsion_id_list_[tor_num].rsd() << " " << torsion_id_list_[tor_num].type() << " " << torsion_id_list_[tor_num].torsion() << " from " << old_tor << " to " << new_tor <<  std::endl;

		// change torsion
		pose.conformation().set_torsion( torsion_id_list_[ tor_num ], new_tor );
	}
}

/// @brief Look at all the bb torsions, make a list of the ones that can move
///
/// I feel like there should be a better way to do this rather than having to hard code so much conectivity information
/// about the residue and the peptide/peptoid checks. It should probably iterate over the move map rather than the pose
/// but done is better than perfect.
void
RandomOmegaFlipMover::setup_torsion_list( core::pose::Pose & pose )
{
	using namespace core;
	using namespace id;

	// clear existing list
	torsion_id_list_.clear();

	// make list
	for ( Size i( 1 ); i <= pose.total_residue(); ++i ) {

		// check to see if peptoid or protein
		if ( pose.residue( i ).type().is_protein() || pose.residue( i ).type().is_peptoid() ) {

			// get omg backbone torsions
			TorsionID omg_tor_id( i, BB, omega_torsion );

			// add moveable torsions to torsion id list
			if ( move_map_->get( omg_tor_id ) ) { torsion_id_list_.push_back( omg_tor_id ); }
		}
	}

	// DEBUG
	//for ( Size i(1); i <= torsion_id_list_.size(); ++i ) {
	//	TR << "DEBUG setup_torsion_list(): " << torsion_id_list_[i].rsd() << " " << torsion_id_list_[i].type() << " " << torsion_id_list_[i].torsion() << std::endl;
	//}
}

protocols::moves::MoverOP
RandomOmegaFlipMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::RandomOmegaFlipMover( *this ) );
}

protocols::moves::MoverOP
RandomOmegaFlipMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::RandomOmegaFlipMover() );
}

void
RandomOmegaFlipMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	if ( !move_map_ ) move_map_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
	protocols::rosetta_scripts::parse_movemap( tag, pose, move_map_, data, false );
}

/// @brief RandomOmegaFlipMoverCreator interface, name of the mover
std::string RandomOmegaFlipMoverCreator::mover_name() {
  return "RandomOmegaFlipMover";
}

/// @brief RandomOmegaFlipMoverCreator interface, returns a unique key name to be used in xml file
std::string RandomOmegaFlipMoverCreator::keyname() const {
  return RandomOmegaFlipMoverCreator::mover_name();
}

/// @brief RandomOmegaFlipMoverCreator interface, return a new instance
protocols::moves::MoverOP RandomOmegaFlipMoverCreator::create_mover() const {
  return protocols::moves::MoverOP( new RandomOmegaFlipMover() );
}


} // simple_moves
} // protocols

/*

TR << "DEBUG" << std::endl;
	move_map_->show( pose.total_residue() );

	TR << "DEBUG TORSION TYPE" << std::endl;
	for( MoveMap::TorsionTypeMap::const_iterator i( move_map_->torsion_type_begin() ), i_end( move_map_->torsion_type_end() ); i != i_end; ++i ) {
		TR << "TORSION TYPE: " << i->first << " " << i->second << std::endl;
	}

	TR << "DEBUG MOVEMAP TORSION ID" << std::endl;
	for( MoveMap::MoveMapTorsionID_Map::const_iterator i( move_map_->movemap_torsion_id_begin() ), i_end( move_map_->movemap_torsion_id_end() ); i != i_end; ++i ) {
		TR << "MM TORSION TYPE: "
		<< i->first.first << " "
		<< i->first.second << " "
		<< i->second << std::endl;
	}

	TR << "DEBUG TORSION ID" << std::endl;
	for( MoveMap::TorsionID_Map::const_iterator i( move_map_->torsion_id_begin() ), i_end( move_map_->torsion_id_end() ); i != i_end; ++i ) {
		TR << "TORSION ID: " << i->first << " " << i->second << std::endl;
	}



 */
