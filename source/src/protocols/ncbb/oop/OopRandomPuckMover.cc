// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ncbb/oop/OopRandomPuckMover.cc
/// @brief OopRandomPuckMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/ncbb/oop/OopRandomPuckMover.fwd.hh>
#include <protocols/ncbb/oop/OopRandomPuckMover.hh>
#include <protocols/ncbb/oop/OopPuckMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
// Package Headers

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
// Random number generator
#include <numeric/random/random.hh>
// Utility Headers
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <core/types.hh>
#include <utility>

// C++ Headers

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.simple_moves.oop.OopRandomPuckMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::chiral;

namespace protocols {
namespace simple_moves {
namespace oop {

/// @details
void OopRandomPuckMover::apply( core::pose::Pose & pose ){

	TR<< "in OopRandomPuckMover::apply" << std::endl;
	//kdrew: for all positions in oop_seq_positions_, input assertion check
	for ( Size i = 1; i <= oop_seq_positions_.size(); i++ ) {
		Size oop_pre_pos = oop_seq_positions_[i];
		Size oop_post_pos = oop_pre_pos+1;
		TR<< "oop_pre_pos:" << oop_pre_pos << " oop_post_pos:" << oop_post_pos << std::endl;

		runtime_assert ( pose.residue(oop_pre_pos).has_variant_type(chemical::OOP_PRE) == 1) ;
		runtime_assert ( pose.residue(oop_post_pos).has_variant_type(chemical::OOP_POST) == 1) ;

		//kdrew: an oop pre position cannot be last position
		runtime_assert ( oop_pre_pos != pose.size() );
		//kdrew: an oop post position cannot be first position
		runtime_assert ( oop_post_pos != 1 );

	}//for

	//kdrew: randomly choose position from oop_seq_positions
	core::Size random_pos = oop_seq_positions_[int(numeric::random::rg().uniform()*oop_seq_positions_.size())+1];

	//kdrew: randomly choose conformation up, down or small angle move
	std::string random_pucker = available_moves_[int(numeric::random::rg().uniform()*available_moves_.size())+1];

	runtime_assert ( random_pucker == "OOP_PUCK_PLUS" || random_pucker == "OOP_PUCK_MINUS" );
	TR << random_pucker <<std::endl;

	oop::OopMoverOP oop_mover;
	ResidueType restype = pose.residue_type( random_pos );

	//kdrew: determine which mover should be used, use D puck movers for chiral D oops
	if ( random_pucker == "OOP_PUCK_PLUS" ) {
		if ( restype.is_d_aa() ) {
			oop_mover = oop::OopMoverOP( new oop::OopDPuckPlusMover( random_pos ) );
		} else {
			oop_mover = oop::OopMoverOP( new oop::OopPuckPlusMover( random_pos ) );
		}
	} else if ( random_pucker == "OOP_PUCK_MINUS" ) {
		if ( restype.is_d_aa() ) {
			oop_mover = oop::OopMoverOP( new oop::OopDPuckMinusMover( random_pos ) );
		} else {
			oop_mover = oop::OopMoverOP( new oop::OopPuckMinusMover( random_pos ) );
		}
	}

	oop_mover->apply( pose );


}//apply

std::string
OopRandomPuckMover::get_name() const {
	return "OopRandomPuckMover";
}

/// @brief
OopRandomPuckMover::OopRandomPuckMover(
) : Mover()
{
	Mover::type( "OopRandomPuckMover" );
}

OopRandomPuckMover::OopRandomPuckMover(
	utility::vector1< core::Size > const & oop_seq_positions
): Mover(), oop_seq_positions_(oop_seq_positions)
{
	Mover::type( "OopRandomPuckMover" );

	available_moves_.push_back("OOP_PUCK_PLUS");
	available_moves_.push_back("OOP_PUCK_MINUS");
}

OopRandomPuckMover::~OopRandomPuckMover()= default;

}//oop
}//simple_moves
}//protocols

