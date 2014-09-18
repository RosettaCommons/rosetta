// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/hbs/HbsRandomSmallMover.cc
/// @brief HbsRandomSmallMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/hbs/HbsRandomSmallMover.fwd.hh>
#include <protocols/simple_moves/hbs/HbsRandomSmallMover.hh>
#include <protocols/simple_moves/hbs/HbsMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
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

/* 
 TODO:
 
 rename to HbsSmallMover, and make the codebase okay with that
 because this isn't actually random, since we have only one HBS position per pose.
 UNLESS we have a vector of HBS positions in the macrocycle... but this seems designed for different OOPs (since you can put an OOP ring anywhere)
 rather than for "randomly choose either OOP_pre_pos or OOP_post_pos to fiddle with
 so the question of how to move different positions on the HBS helix differently
 seems to have a solution not found here.
 
*/

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.simple_moves.hbs.HbsRandomSmallMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace hbs {

///@details
void HbsRandomSmallMover::apply( core::pose::Pose & pose ){

	using numeric::conversions::radians;
	using numeric::conversions::degrees;
	
	TR<< "in HbsRandomSmallMover::apply" << std::endl;
	//kdrew: input assertion check
	Size hbs_pre_pos = hbs_seq_position_;
	Size hbs_post_pos = hbs_pre_pos+2;                   
	TR<< "hbs_pre_pos:" << hbs_pre_pos << " hbs_post_pos:" << hbs_post_pos << std::endl;

	//runtime_assert ( pose.residue(hbs_pre_pos).has_variant_type(chemical::HBS_PRE) == 1) ;
	runtime_assert ( pose.residue(hbs_post_pos).has_variant_type(chemical::HBS_POST) == 1) ;
	//kdrew: an hbs pre position cannot be last position
	runtime_assert ( hbs_pre_pos != pose.total_residue() );
	//kdrew: an hbs post position cannot be first position
	runtime_assert ( hbs_post_pos != 1 );

	//kdrew: randomly choose position from hbs_seq_positions
	//core::Size random_pos = hbs_post_pos + int(numeric::random::rg().uniform()*(hbs_length_-3))+1;
	core::Size random_pos = hbs_pre_pos+int(numeric::random::rg().uniform()*(hbs_length_-1))+1;
	
	hbs::HbsMoverOP hbs_mover ( new hbs::HbsMover( random_pos ) );
	Real small_angle = max_small_angle_/2.0; ///< this is max_angle/2, which is the deviation from the angle input
	Real phi_angle = basic::periodic_range( pose.phi( random_pos ) - small_angle + numeric::random::rg().uniform() * max_small_angle_, 360.0 );
	//kdrew: no phi angle for n-terms, angle that gets changed is CYH-N-Ca-C
	if( pose.residue_type( random_pos ).is_lower_terminus() )
	{ 
		AtomID aidCYH( pose.residue(random_pos).atom_index("CYH"), random_pos );
		AtomID aidN( pose.residue(random_pos).atom_index("N"), random_pos );
		AtomID aidCA( pose.residue(random_pos).atom_index("CA"), random_pos );
		AtomID aidC( pose.residue(random_pos).atom_index("C"), random_pos );

		Real CYH_N_Ca_C_angle = degrees( pose.conformation().torsion_angle( aidCYH, aidN, aidCA, aidC ) ); 
        phi_angle = basic::periodic_range( CYH_N_Ca_C_angle - small_angle + numeric::random::rg().uniform() * max_small_angle_, 360.0 ) - 180.0; 
	}

	Real psi_angle = basic::periodic_range( pose.psi( random_pos ) - small_angle + numeric::random::rg().uniform() * max_small_angle_, 360.0 );
	hbs_mover->set_phi( phi_angle );
	hbs_mover->set_psi( psi_angle );
	hbs_mover->apply( pose );
	
}//apply

std::string
HbsRandomSmallMover::get_name() const {
	return "HbsRandomSmallMover";
}

///@brief
/*
 awatkins: fuck the empty constructor. why would we use it? why is it throwing an error?
 
 HbsRandomSmallMover::HbsRandomSmallMover(
) : Mover()
{
	Mover::type( "HbsRandomSmallMover" );
}*/

HbsRandomSmallMover::HbsRandomSmallMover( core::Size hbs_seq_position, core::Real max_small_angle): Mover(), hbs_seq_position_(hbs_seq_position), hbs_length_(12), max_small_angle_(max_small_angle)
{
	Mover::type( "HbsRandomSmallMover" );
}

HbsRandomSmallMover::HbsRandomSmallMover( core::Size hbs_seq_position/*, core::Size hbs_length*/ ): Mover(), hbs_seq_position_(hbs_seq_position)/*, hbs_length_(hbs_length)*/, max_small_angle_(2.0)
{
	Mover::type( "HbsRandomSmallMover" );
}

/*HbsRandomSmallMover::HbsRandomSmallMover( core::Size  hbs_seq_position, core::Size hbs_length, core::Real max_small_angle ): Mover(), hbs_length_(hbs_length), hbs_seq_position_(hbs_seq_position), max_small_angle_( max_small_angle )
{
	Mover::type( "HbsRandomSmallMover" );
}*/

HbsRandomSmallMover::~HbsRandomSmallMover(){}


}//hbs
}//simple_moves
}//protocols

