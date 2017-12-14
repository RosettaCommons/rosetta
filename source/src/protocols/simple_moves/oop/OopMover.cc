// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/oop/OopMover.cc
/// @brief OopMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/oop/OopMover.hh>
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

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.simple_moves.oop.OopMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace oop {

/*
kdrew: the helper function moves a single oop with the given move
pose: pose to make change to
oop_pre_pos: position of first residue of oop to make change to
phi_angle: value to change phi angle to
psi_angle: value to change psi angle to
*/
void OopMover::apply( core::pose::Pose & pose )
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	runtime_assert ( pose.residue(oop_pre_pos_).has_variant_type(chemical::OOP_PRE) == 1) ;
	runtime_assert ( pose.residue(oop_post_pos_).has_variant_type(chemical::OOP_POST) == 1) ;
	//kdrew: an oop pre position cannot be last position
	runtime_assert ( oop_pre_pos_ != pose.size() );
	//kdrew: an oop post position cannot be first position
	runtime_assert ( oop_post_pos_ != 1 );

	//kdrew: first residue is a special case, no phi torsion
	if ( pose.residue_type( oop_pre_pos_ ).is_lower_terminus() ) {
		//kdrew: no phi angle, adjust CYP-N-Ca-C torsion
		Real offset = 180 + phi_angle_ ;

		AtomID aidCYP( pose.residue(oop_pre_pos_).atom_index("CYP"), oop_pre_pos_ );
		AtomID aidN( pose.residue(oop_pre_pos_).atom_index("N"), oop_pre_pos_ );
		AtomID aidCA( pose.residue(oop_pre_pos_).atom_index("CA"), oop_pre_pos_ );
		AtomID aidC( pose.residue(oop_pre_pos_).atom_index("C"), oop_pre_pos_ );

		pose.conformation().set_torsion_angle( aidCYP, aidN, aidCA, aidC, radians( offset ) );
	}


	TR << "current oop_pre ("<< oop_pre_pos_<< ") PHI: " << pose.phi(oop_pre_pos_) << std::endl;
	pose.set_phi(oop_pre_pos_, phi_angle_);
	TR << "new oop_pre ("<< oop_pre_pos_<< ") PHI: " << pose.phi(oop_pre_pos_) << std::endl;
	//pose.dump_pdb( "rosetta_out_oop_pre.pdb" );
	TR << "current oop_pre ("<< oop_pre_pos_<< ") PSI: " << pose.psi(oop_pre_pos_) << std::endl;
	pose.set_psi(oop_pre_pos_, psi_angle_);
	TR << "new oop_pre ("<< oop_pre_pos_<< ") PSI: " << pose.psi(oop_pre_pos_) << std::endl;
	//pose.dump_pdb( "rosetta_out_oop_pre_pos_t.pdb" );

	update_hydrogens_( pose );

}

std::string
OopMover::get_name() const {
	return "OopMover";
}

/// @brief
OopMover::OopMover(
	core::Size oop_seq_position
): Mover(), oop_pre_pos_(oop_seq_position), oop_post_pos_(oop_seq_position+1)
{
	Mover::type( "OopMover" );
}

OopMover::OopMover(
	core::Size oop_seq_position,
	core::Real phi_angle,
	core::Real psi_angle
): Mover(), oop_pre_pos_(oop_seq_position), oop_post_pos_(oop_seq_position+1), phi_angle_(phi_angle), psi_angle_(psi_angle)
{
	Mover::type( "OopMover" );

}

OopMover::~OopMover()= default;

void OopMover::update_hydrogens_( core::pose::Pose & pose )
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	AtomID aidVZP( pose.residue(oop_pre_pos_).atom_index("VZP"), oop_pre_pos_);
	AtomID aidCYP( pose.residue(oop_pre_pos_).atom_index("CYP"), oop_pre_pos_);
	AtomID aidN( pose.residue(oop_pre_pos_).atom_index("N"), oop_pre_pos_);
	AtomID aidCA( pose.residue(oop_pre_pos_).atom_index("CA"), oop_pre_pos_);
	AtomID aid1HYP( pose.residue(oop_pre_pos_).atom_index("1HYP"), oop_pre_pos_);
	AtomID aidVYP( pose.residue(oop_post_pos_).atom_index("VYP"), oop_post_pos_);
	AtomID aidCZP( pose.residue(oop_post_pos_).atom_index("CZP"), oop_post_pos_);
	AtomID aidN_Z( pose.residue(oop_post_pos_).atom_index("N"), oop_post_pos_);
	AtomID aidCA_Z( pose.residue(oop_post_pos_).atom_index("CA"), oop_post_pos_);
	AtomID aid1HZP( pose.residue(oop_post_pos_).atom_index("1HZP"), oop_post_pos_);

	//kdrew: definitions for xyz coordinates
	Vector const& CA_xyz ( pose.residue(oop_pre_pos_).xyz("CA") );
	Vector const& N_xyz ( pose.residue(oop_pre_pos_).xyz("N") );
	Vector const& CYP_xyz ( pose.residue(oop_pre_pos_).xyz("CYP") );
	Vector const& CAZ_xyz ( pose.residue(oop_post_pos_).xyz("CA") );
	Vector const& NZ_xyz ( pose.residue(oop_post_pos_).xyz("N") );
	Vector const& CZP_xyz ( pose.residue(oop_post_pos_).xyz("CZP") );


	//kdrew: VZP and VYP are relative to the 1Hs, so move the 1Hs and VZP|VYP and 2H moves with it
	//kdrew: calculate the difference between the virtual atom VZP torsion and the torsion involving the real CZP atom (on the post residue)
	Real VZP_torsion_correction = numeric::dihedral_degrees( CZP_xyz, CYP_xyz, N_xyz, CA_xyz ) - degrees( pose.conformation().torsion_angle( aidVZP, aidCYP, aidN, aidCA ) );
	//kdrew: change the 1HYP torsion by the corrected amount
	Real torsion_1HYP = degrees(pose.conformation().torsion_angle(aid1HYP,aidCYP,aidN,aidCA)) + VZP_torsion_correction;
	pose.conformation().set_torsion_angle(aid1HYP,aidCYP,aidN,aidCA,radians(torsion_1HYP));
	//pose.dump_pdb( "rosetta_out_oop_pre_mvH.pdb" );

	//kdrew: calculate the difference between the virtual atom VYP torsion and the torsion involving the real CYP atom (on the pre residue)
	Real VYP_torsion_correction = numeric::dihedral_degrees( CYP_xyz, CZP_xyz, NZ_xyz, CAZ_xyz ) - degrees( pose.conformation().torsion_angle( aidVYP, aidCZP, aidN_Z, aidCA_Z ));
	//kdrew: change the 1HZP torsion by the corrected amount
	Real torsion_1HZP = degrees(pose.conformation().torsion_angle(aid1HZP,aidCZP,aidN_Z,aidCA_Z)) + VYP_torsion_correction;
	pose.conformation().set_torsion_angle(aid1HZP,aidCZP,aidN_Z,aidCA_Z,radians(torsion_1HZP));
	//pose.dump_pdb( "rosetta_out_oop_pre_pos_t_mvH.pdb" );

	//kdrew: possibly minimize pose here

	//kdrew: debug printing
	//kdrew: VZP CYP N CA
	TR<<"VZP-CYP-N-CA: "<< degrees(pose.conformation().torsion_angle(aidVZP,aidCYP,aidN,aidCA)) <<std::endl;
	//kdrew: CZP CYP N CA
	TR<<"CZP-CYP-N-CA: "<< numeric::dihedral_degrees(CZP_xyz,CYP_xyz,N_xyz,CA_xyz) <<std::endl;
	//kdrew: 1HYP CYP N CA
	TR<<"1HYP-CYP-N-CA: "<< degrees(pose.conformation().torsion_angle(aid1HYP,aidCYP,aidN,aidCA)) <<std::endl;

	//kdrew: VYP CZP N_Z CA_Z
	TR<<"VYP-CZP-N-CA: "<< degrees(pose.conformation().torsion_angle(aidVYP,aidCZP,aidN_Z,aidCA_Z)) << std::endl;
	//kdrew: CYP CZP N_Z CA_Z
	TR<<"CYP-CZP-N-CA: "<< numeric::dihedral_degrees(CYP_xyz,CZP_xyz,NZ_xyz,CAZ_xyz) <<std::endl;
	//kdrew: 1HZP CZP N_Z CA_Z
	TR<<"1HZP-CZP-N-CA: "<< degrees(pose.conformation().torsion_angle(aid1HZP,aidCZP,aidN_Z,aidCA_Z)) << std::endl;

}

}//oop
}//simple_moves
}//protocols

