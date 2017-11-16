// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/residue_selector/util.hh>

// utility headers
#include <utility>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// numeric headers
#include <numeric/random/random.hh>

// basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.simple_moves.RandomTorsionMover" );

namespace protocols {
namespace simple_moves {


RandomTorsionMover::RandomTorsionMover() :
	Mover("RandomTorsionMover"),
	movemap_factory_( /* nullptr */ ),
	move_map_( /* nullptr */ ),
	max_angle_( 0 ),
	num_moves_( 0 )
{}

RandomTorsionMover::RandomTorsionMover( core::select::movemap::MoveMapFactoryCOP movemap_factory, core::Real max_angle, core::Size num_moves ) :
	Mover("RandomTorsionMover"),
	movemap_factory_(std::move( movemap_factory )),
	max_angle_( max_angle ),
	num_moves_( num_moves )
{}

RandomTorsionMover::RandomTorsionMover( core::kinematics::MoveMapOP move_map, core::Real max_angle, core::Size num_moves ) :
	Mover("RandomTorsionMover"),
	move_map_(std::move( move_map )),
	max_angle_( max_angle ),
	num_moves_( num_moves )
{}

RandomTorsionMover::RandomTorsionMover( RandomTorsionMover const & other ) :
	Mover("RandomTorsionMover"),
	movemap_factory_( other.movemap_factory_),
	max_angle_( other.max_angle_ ),
	num_moves_( other.num_moves_ )
{
	if ( other.move_map_ ) { // Don't clone null ptr
		move_map_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap( *other.move_map_ ) );
	}
}


RandomTorsionMover::~RandomTorsionMover()= default;

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
			Size tor_num( numeric::random::rg().random_range( 1, torsion_id_list_.size() ) );

			// calc randomly purturbed value
			Real old_tor( pose.conformation().torsion( torsion_id_list_[ tor_num ] ) );
			Real new_tor( periodic_range( old_tor - ( max_angle_ / 2 ) + ( numeric::random::rg().uniform() * max_angle_ ), 360.0 ) );

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

	core::kinematics::MoveMapCOP move_map( movemap( pose ) );
	// make list
	for ( Size i( 1 ); i <= pose.size(); ++i ) {

		// check to see if peptoid or protein
		if ( pose.residue( i ).type().is_protein() || pose.residue( i ).type().is_peptoid() ) {

			// get three backbone torsions
			TorsionID phi_tor_id( i, BB, phi_torsion );
			TorsionID psi_tor_id( i, BB, psi_torsion );
			TorsionID omg_tor_id( i, BB, omega_torsion );

			// add moveable torsions to torsion id list
			if ( move_map->get( phi_tor_id ) ) { torsion_id_list_.push_back( phi_tor_id ); }
			if ( move_map->get( psi_tor_id ) ) { torsion_id_list_.push_back( psi_tor_id ); }
			if ( move_map->get( omg_tor_id ) ) { torsion_id_list_.push_back( omg_tor_id ); }
		}
	}

	// DEBUG
	//for ( Size i(1); i <= torsion_id_list_.size(); ++i ) {
	// TR << "DEBUG torsion list" << torsion_id_list_[i].rsd() << " " << torsion_id_list_[i].type() << " " << torsion_id_list_[i].torsion() << std::endl;
	//}
}

protocols::moves::MoverOP
RandomTorsionMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::RandomTorsionMover( *this ) );
}

protocols::moves::MoverOP
RandomTorsionMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::RandomTorsionMover() );
}

void
RandomTorsionMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	max_angle_ = tag->getOption< core::Real >( "max_angle", max_angle_ );
	num_moves_ = tag->getOption< core::Size >( "num_moves", num_moves_ );

	movemap_factory( protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data ) );
}

/// @brief RandomTorsionMoverCreator interface, name of the mover
// XRW TEMP std::string RandomTorsionMover::mover_name() {
// XRW TEMP  return "RandomTorsionMover";
// XRW TEMP }

/// @brief RandomTorsionMoverCreator interface, returns a unique key name to be used in xml file
// XRW TEMP std::string RandomTorsionMoverCreator::keyname() const {
// XRW TEMP  return RandomTorsionMover::mover_name();
// XRW TEMP }

/// @brief RandomTorsionMoverCreator interface, return a new instance
// XRW TEMP protocols::moves::MoverOP RandomTorsionMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RandomTorsionMover() );
// XRW TEMP }

std::string RandomTorsionMover::get_name() const {
	return mover_name();
}

std::string RandomTorsionMover::mover_name() {
	return "RandomTorsionMover";
}

core::kinematics::MoveMapCOP
RandomTorsionMover::movemap( core::pose::Pose const & pose ) {
	if ( move_map_ ) {
		return move_map_;
	} else if ( movemap_factory_ ) {
		return movemap_factory_->create_movemap_from_pose( pose );
	} else {
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		return movemap;
	}
}

void RandomTorsionMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "max_angle", xsct_real, "XRW_TODO" )
		+ XMLSchemaAttribute( "num_moves", xsct_non_negative_integer, "XRW_TODO" );

	XMLSchemaSimpleSubelementList subelements;
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, subelements );
	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW_TODO", attlist, subelements );
}

std::string RandomTorsionMoverCreator::keyname() const {
	return RandomTorsionMover::mover_name();
}

protocols::moves::MoverOP
RandomTorsionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RandomTorsionMover );
}

void RandomTorsionMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RandomTorsionMover::provide_xml_schema( xsd );
}


} // simple_moves
} // protocols

/*

TR << "DEBUG" << std::endl;
move_map_->show( pose.size() );

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
