// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/triazolamer/TriazolamerRandomSmallMover.cc
/// @brief TriazolamerRandomSmallMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/triazolamer/TriazolamerRandomSmallMover.fwd.hh>
#include <protocols/simple_moves/triazolamer/TriazolamerRandomSmallMover.hh>
#include <protocols/simple_moves/triazolamer/TriazolamerMover.hh>
// Package Headers

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
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

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.triazolamer.TriazolamerRandomSmallMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace triazolamer {

/// @details
void TriazolamerRandomSmallMover::apply( core::pose::Pose & pose ){

	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	TR<< "in TriazolamerRandomSmallMover::apply" << std::endl;
	//kdrew: for all positions in triazolamer_seq_positions_, input assertion check
	for ( Size i = 1; i <= triazolamer_seq_positions_.size(); i++ ) {
		Size triazolamer_pre_pos = triazolamer_seq_positions_[i];
		Size triazolamer_post_pos = triazolamer_pre_pos+1;
		TR<< "triazolamer_pre_pos:" << triazolamer_pre_pos << " triazolamer_post_pos:" << triazolamer_post_pos << std::endl;

		runtime_assert ( pose.residue(triazolamer_pre_pos).has_variant_type(chemical::TRIAZOLAMERC) == 1) ;
		runtime_assert ( pose.residue(triazolamer_post_pos).has_variant_type(chemical::TRIAZOLAMERN) == 1) ;
		//kdrew: an triazolamer pre position cannot be last position
		runtime_assert ( triazolamer_pre_pos != pose.total_residue() );
		//kdrew: an triazolamer post position cannot be first position
		runtime_assert ( triazolamer_post_pos != 1 );

	}//for


	//kdrew: randomly choose position from triazolamer_seq_positions
	core::Size random_pos = triazolamer_seq_positions_[int(numeric::random::rg().uniform()*triazolamer_seq_positions_.size())+1];

	triazolamer::TriazolamerMoverOP triazolamer_mover( new triazolamer::TriazolamerMover( random_pos ) );
	Real small_angle = max_small_angle_/2.0; ///< this is max_angle/2, which is the deviation from the angle input
	Real phi_angle = basic::periodic_range( pose.phi( random_pos ) - small_angle + numeric::random::rg().uniform() * max_small_angle_, 360.0 );
	//kdrew: no phi angle for n-terms, angle that gets changed is CYP-N-Ca-C
	if ( pose.residue_type( random_pos ).is_lower_terminus() ) {
		AtomID aidCYP( pose.residue(random_pos).atom_index("CYP"), random_pos );
		AtomID aidN( pose.residue(random_pos).atom_index("N"), random_pos );
		AtomID aidCA( pose.residue(random_pos).atom_index("CA"), random_pos );
		AtomID aidC( pose.residue(random_pos).atom_index("C"), random_pos );

		Real CYP_N_Ca_C_angle = degrees( pose.conformation().torsion_angle( aidCYP, aidN, aidCA, aidC ) );
		phi_angle = basic::periodic_range( CYP_N_Ca_C_angle - small_angle + numeric::random::rg().uniform() * max_small_angle_, 360.0 ) - 180.0;
	}

	Real psi_angle = basic::periodic_range( pose.psi( random_pos ) - small_angle + numeric::random::rg().uniform() * max_small_angle_, 360.0 );
	triazolamer_mover->set_phi( phi_angle );
	triazolamer_mover->set_psi( psi_angle );
	triazolamer_mover->apply( pose );
	//triazolamer_puck_mover_helper( pose, random_pos, phi_angle, psi_angle );

}//apply

std::string
TriazolamerRandomSmallMover::get_name() const {
	return "TriazolamerRandomSmallMover";
}

/// @brief
TriazolamerRandomSmallMover::TriazolamerRandomSmallMover(
) : Mover()
{
	Mover::type( "TriazolamerRandomSmallMover" );
}

TriazolamerRandomSmallMover::TriazolamerRandomSmallMover( utility::vector1< core::Size > triazolamer_seq_positions ): Mover(), triazolamer_seq_positions_(triazolamer_seq_positions), max_small_angle_(0.0)
{
	Mover::type( "TriazolamerRandomSmallMover" );
}

TriazolamerRandomSmallMover::TriazolamerRandomSmallMover( utility::vector1< core::Size > triazolamer_seq_positions, core::Real max_small_angle ): Mover(), triazolamer_seq_positions_(triazolamer_seq_positions), max_small_angle_( max_small_angle )
{
	Mover::type( "TriazolamerRandomSmallMover" );
}

TriazolamerRandomSmallMover::~TriazolamerRandomSmallMover(){}


}//triazolamer
}//simple_moves
}//protocols

