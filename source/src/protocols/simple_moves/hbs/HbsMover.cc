// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/hbs/HbsMover.cc
/// @brief HbsMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/hbs/HbsMover.hh>
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
// Utility Headers
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <core/types.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.simple_moves.hbs.HbsMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace hbs {
// awatkins: This should NEVER be called on residues that actually have HBS atom labels, fixing... 7/14/13
/*
kdrew: the helper function moves a single hbs with the given move
pose: pose to make change to
hbs_pre_pos: position of first residue of hbs to make change to
phi_angle: value to change phi angle to
psi_angle: value to change psi angle to
*/
void HbsMover::apply( core::pose::Pose & pose )
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	//runtime_assert ( pose.residue(hbs_pre_pos_).has_variant_type(chemical::HBS_PRE) == 1) ;
	//runtime_assert ( pose.residue(hbs_post_pos_).has_variant_type(chemical::HBS_POST) == 1) ;
	//kdrew: an hbs pre position cannot be last position
	//runtime_assert ( hbs_pre_pos_ != pose.total_residue() );
	//kdrew: an hbs post position cannot be first position
	//runtime_assert ( hbs_post_pos_ != 1 );

	//kdrew: first residue is a special case, no phi torsion
	/*if( pose.residue_type( hbs_pre_pos_ ).is_lower_terminus() )
	{
		//kdrew: no phi angle, adjust CYP-N-Ca-C torsion
		Real offset = 180 + phi_angle_ ;

		AtomID aidCYH( pose.residue(hbs_pre_pos_).atom_index("CYH"), hbs_pre_pos_ );
		AtomID aidN( pose.residue(hbs_pre_pos_).atom_index("N"), hbs_pre_pos_ );
		AtomID aidCA( pose.residue(hbs_pre_pos_).atom_index("CA"), hbs_pre_pos_ );
		AtomID aidC( pose.residue(hbs_pre_pos_).atom_index("C"), hbs_pre_pos_ );

		pose.conformation().set_torsion_angle( aidCYH, aidN, aidCA, aidC, radians( offset ) ); 
	}*/


	TR << "current pos ("<< seq_pos_<< ") PHI: " << pose.phi(seq_pos_) << std::endl;
	pose.set_phi(seq_pos_, phi_angle_);
	TR << "new pos ("<< seq_pos_<< ") PHI: " << pose.phi(seq_pos_) << std::endl;
	//pose.dump_pdb( "rosetta_out_hbs_pre.pdb" );
	TR << "current pos ("<< seq_pos_<< ") PSI: " << pose.psi(seq_pos_) << std::endl;
	pose.set_psi(seq_pos_, psi_angle_);
	TR << "new pos ("<< seq_pos_<< ") PSI: " << pose.psi(seq_pos_) << std::endl;
	//pose.dump_pdb( "rosetta_out_hbs_pre_pos_t.pdb" );

}

std::string
HbsMover::get_name() const {
	return "HbsMover";
}

///@brief
HbsMover::HbsMover( 
		core::Size seq_position 
	): Mover(), seq_pos_(seq_position)
{
	Mover::type( "HbsMover" );
}

HbsMover::HbsMover( 
		core::Size seq_position, 
		core::Real phi_angle,
		core::Real psi_angle
	): Mover(), seq_pos_(seq_position), phi_angle_(phi_angle), psi_angle_(psi_angle)
{
	Mover::type( "HbsMover" );

}

HbsMover::~HbsMover(){}

}//hbs
}//simple_moves
}//protocols

