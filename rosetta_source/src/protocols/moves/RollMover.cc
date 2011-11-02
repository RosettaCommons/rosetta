// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/RollMover.cc
/// @brief RollMover methods implemented
/// @author

// Unit Headers
#include <protocols/moves/RollMover.hh>
#include <protocols/moves/RollMoverCreator.hh>
// Package Headers

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyz.functions.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
// AUTO-REMOVED #include <numeric/xyz.io.hh>
// Random number generator
#include <numeric/random/random.hh>
// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>


// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;
static numeric::random::RandomGenerator RG(456732);
static basic::Tracer TR( "protocols.moves.RollMover" );

using namespace core;

namespace protocols {
namespace moves {

///@details
void RollMover::apply( core::pose::Pose & pose ){

	utility::vector1< utility::vector1< numeric::xyzVector< core::Real> > >  coords; // some of these will change for sure

	// now we have a vector1 of vector1s with a numeric::xyzVector
	// access will look like coords[residue][atom] to get xyz

	core::Size const nres( pose.total_residue() );
  coords.resize( nres );

  for ( Size i=1; i<= nres; ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		core::Size const number_atoms_this_residue( rsd.natoms() );
    if ( number_atoms_this_residue ) {
      coords[i].resize( number_atoms_this_residue );
      for ( Size j=1; j <= number_atoms_this_residue; ++j ) {
        coords[i][j] = rsd.atom( j ).xyz();
			}
		}
	}

	angle_ = min_angle_ + ( max_angle_ - min_angle_ ) * RG.uniform();
	numeric::xyzMatrix< core::Real > rotation_matrix( numeric::rotation_matrix_degrees(axis_, angle_ ) );
	//move to origin
	for ( core::Size i =start_res_; i <= stop_res_; ++i ) {
		for ( core::Size j = 1; j <= coords[i].size(); ++j ) {

			// this may look strange but in a global coordinate system
			// rotation about an axis is easily done by movement to the origin
			// rotation and then movement back

			coords[i][j] = coords[i][j] - translate_; // translate to origin
      coords[i][j] = rotation_matrix * coords[i][j]; // rotate atom
      coords[i][j] = coords[i][j] + translate_; // reverse translate

    }
  }

	// now update pose with new coordinates
  for ( core::Size i =start_res_; i <= stop_res_; ++i ) {
		for ( core::Size j = 1; j <= coords[i].size(); ++j ) {
			id::AtomID id( j, i );
			pose.set_xyz( id, coords[i][j]);
		}
	}



}//apply


void
RollMover::set_min_max_angles( core::Real min_angle, core::Real max_angle ) {
	min_angle_ = min_angle;
	max_angle_ = max_angle;
}

std::string
RollMover::get_name() const {
	return "RollMover";
}

void
RollMover::parse_my_tag( 
	utility::tag::TagPtr const tag,
	DataMap & /*datamap*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/ )
{       

        /*parse start_res*/
	if( tag->hasOption("start_res") ) {
		start_res_ = tag->getOption<core::Size>("start_res");
	} else {
		utility_exit_with_message("RollMover requires start_res option");
	}

	/*parse stop_res*/
	if( tag->hasOption("stop_res") ) {
		stop_res_ = tag->getOption<core::Size>("stop_res");
	} else {
		utility_exit_with_message("RollMover requires stop_res option");
	}

	/*parse min_angle*/
	if( tag->hasOption("min_angle") ) {
		min_angle_ = tag->getOption<core::Real>("min_angle");
	} else {
		utility_exit_with_message("RollMover requires min_angle option");
	}

	/*parse max_angle*/
	if( tag->hasOption("max_angle") ) {
		max_angle_ = tag->getOption<core::Real>("max_angle");
	} else {
		utility_exit_with_message("RollMover requires max_angle option");
	}

	bool axis_option_parsed = false;
	bool translate_option_parsed = false;

	foreach( utility::tag::TagPtr const child_tag, tag->getTags() ){
		std::string name= child_tag->getName();

		if( name == "axis" ) {
			/*parse axis x,y,z*/
			axis_ = protocols::rosetta_scripts::parse_xyz_vector(child_tag);
			axis_option_parsed = true;
		} else if ( name == "translate") {
			/*parse translate x,y,z*/
			translate_ = protocols::rosetta_scripts::parse_xyz_vector(child_tag);
			translate_option_parsed = true;
		}

	}

	if ( !axis_option_parsed ) {
		utility_exit_with_message("RollMover requires axis option");
	}
	if ( !translate_option_parsed ) {
		utility_exit_with_message("RollMover requires translate option");
	}
	
}

std::string
RollMoverCreator::keyname() const
{
	return RollMoverCreator::mover_name();
}

std::string
RollMoverCreator::mover_name()
{
	return "RollMover";
}


protocols::moves::MoverOP
RollMoverCreator::create_mover() const {
	return new RollMover;
}

///@brief required in the context of the parser/scripting scheme
MoverOP
RollMover::fresh_instance() const
{
	return new RollMover;
}

///@brief required in the context of the parser/scripting scheme
MoverOP
RollMover::clone() const
{
	return new RollMover( *this );
}

///@brief
RollMover::RollMover(
) : Mover()
{
	Mover::type( "RollMover" );
}

	RollMover::RollMover( core::Size start_res, core::Size stop_res, core::Real min_angle, core::Real max_angle, numeric::xyzVector< core::Real > axis, numeric::xyzVector< core::Real > translate
												): Mover(), start_res_(start_res), stop_res_(stop_res), min_angle_(min_angle), max_angle_(max_angle), axis_(axis), translate_(translate)
{
	Mover::type( "RollMover" );
}

RollMover::~RollMover(){}

}//moves
}//protocols

