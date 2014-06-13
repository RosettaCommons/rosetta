// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/RandomTorsionMover.cc
/// @brief RandomTorsionMover methods implemented
/// @author

// unit headers
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/RandomTorsionMoverCreator.hh>

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

static basic::Tracer TR( "protocols.simple_moves.RandomTorsionMover" );

namespace protocols {
namespace simple_moves {

static numeric::random::RandomGenerator RG(12345678);

RandomTorsionMover::RandomTorsionMover() :
	Mover("RandomTorsionMover"),
	move_map_( 0 ),
	max_angle_( 0 ),
	num_moves_( 0 )
{}

RandomTorsionMover::RandomTorsionMover( core::kinematics::MoveMapOP move_map, core::Real max_angle, core::Size num_moves = 1) :
	Mover("RandomTorsionMover"),
	move_map_( move_map ),
	max_angle_( max_angle ),
	num_moves_( num_moves )
{}

RandomTorsionMover::RandomTorsionMover( RandomTorsionMover const & other ) :
	Mover("RandomTorsionMover"),
	move_map_( new core::kinematics::MoveMap( *other.move_map_ ) ),
	max_angle_( other.max_angle_ ),
	num_moves_( other.num_moves_ )
{}


RandomTorsionMover::~RandomTorsionMover(){}

void
RandomTorsionMover::apply( core::pose::Pose & pose )
{
	using namespace basic;
	using namespace core;
	using namespace conformation;
	using namespace kinematics;

	// setup the torsion id list
	setup_torsion_list( pose );

	// empty move maps have zero torsions
	if ( torsion_id_list_.size() ) {

		for ( Size i(1); i <= num_moves_; ++i ) {

			// randomly select a free torsion
			Size tor_num( RG.random_range( 1, torsion_id_list_.size() ) );

			// calc randomly purturbed value
			Real old_tor( pose.conformation().torsion( torsion_id_list_[ tor_num ] ) );
			Real new_tor( periodic_range( old_tor - ( max_angle_ / 2 ) + ( RG.uniform() * max_angle_ ), 360.0 ) );

			//TR << "DEBUG: Setting torsion " << torsion_id_list_[tor_num].rsd() << " " << torsion_id_list_[tor_num].type() << " " << torsion_id_list_[tor_num].torsion() << " from " << old_tor << " to " << new_tor <<  std::endl;

			// change torsion
			pose.conformation().set_torsion( torsion_id_list_[ tor_num ], new_tor );
		}
	}
}

/// @brief Look at all the bb torsions, make a list of the ones that can move
/// I feel like there should be a better way to do this rather than having to hard code so much conectivity information
/// about the residue and the peptide/peptoid checks.
void
RandomTorsionMover::setup_torsion_list( core::pose::Pose & pose )
{
	using namespace core;
	using namespace id;

	// clear existing list
	torsion_id_list_.clear();

	// make list
	for ( Size i( 1 ); i <= pose.total_residue(); ++i ) {

		// check to see if peptoid or protein
		if ( pose.residue( i ).type().is_protein() || pose.residue( i ).type().is_peptoid() ) {

			// get three backbone torsions
			TorsionID phi_tor_id( i, BB, phi_torsion );
			TorsionID psi_tor_id( i, BB, psi_torsion );
			TorsionID omg_tor_id( i, BB, omega_torsion );

			// add moveable torsions to torsion id list
			if ( move_map_->get( phi_tor_id ) ) { torsion_id_list_.push_back( phi_tor_id ); }
			if ( move_map_->get( psi_tor_id ) ) { torsion_id_list_.push_back( psi_tor_id ); }
			if ( move_map_->get( omg_tor_id ) ) { torsion_id_list_.push_back( omg_tor_id ); }
		}
	}

	// DEBUG
	//for ( Size i(1); i <= torsion_id_list_.size(); ++i ) {
	//	TR << "DEBUG torsion list" << torsion_id_list_[i].rsd() << " " << torsion_id_list_[i].type() << " " << torsion_id_list_[i].torsion() << std::endl;
	//}
}

protocols::moves::MoverOP
RandomTorsionMover::clone() const
{
	return new protocols::simple_moves::RandomTorsionMover( *this );
}

protocols::moves::MoverOP
RandomTorsionMover::fresh_instance() const
{
	return new protocols::simple_moves::RandomTorsionMover();
}

void
RandomTorsionMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	if ( !move_map_ ) move_map_ = new core::kinematics::MoveMap();

	max_angle_ = tag->getOption< core::Real >( "max_angle", max_angle_ );
	num_moves_ = tag->getOption< core::Size >( "num_moves", num_moves_ );

	protocols::rosetta_scripts::parse_movemap( tag, pose, move_map_, data, false );
}

/// @brief RandomTorsionMoverCreator interface, name of the mover
std::string RandomTorsionMoverCreator::mover_name() {
  return "RandomTorsionMover";
}

/// @brief RandomTorsionMoverCreator interface, returns a unique key name to be used in xml file
std::string RandomTorsionMoverCreator::keyname() const {
  return RandomTorsionMoverCreator::mover_name();
}

/// @brief RandomTorsionMoverCreator interface, return a new instance
protocols::moves::MoverOP RandomTorsionMoverCreator::create_mover() const {
  return new RandomTorsionMover();
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
