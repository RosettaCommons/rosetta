// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/triazolamer/TriazolamerMover.cc
/// @brief TriazolamerMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
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
// Utility Headers
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <core/types.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.simple_moves.triazolamer.TriazolamerMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace triazolamer {

/*
kdrew: the helper function moves a single triazolamer with the given move
pose: pose to make change to
triazolamer_pre_pos: position of first residue of triazolamer to make change to
phi_angle: value to change phi angle to
psi_angle: value to change psi angle to
*/
void TriazolamerMover::apply( core::pose::Pose & pose )
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	runtime_assert ( pose.residue(triazolamer_pre_pos_).has_variant_type(chemical::TRIAZOLAMERC) == 1) ;
	runtime_assert ( pose.residue(triazolamer_post_pos_).has_variant_type(chemical::TRIAZOLAMERN) == 1) ;
	//kdrew: an triazolamer pre position cannot be last position
	runtime_assert ( triazolamer_pre_pos_ != pose.total_residue() );
	//kdrew: an triazolamer post position cannot be first position
	runtime_assert ( triazolamer_post_pos_ != 1 );

	//kdrew: first residue is a special case, no phi torsion
	if ( pose.residue_type( triazolamer_pre_pos_ ).is_lower_terminus() ) {
		//kdrew: no phi angle, adjust CYP-N-Ca-C torsion
		Real offset = 180 + phi_angle_ ;

		AtomID aidCYP( pose.residue(triazolamer_pre_pos_).atom_index("CYP"), triazolamer_pre_pos_ );
		AtomID aidN( pose.residue(triazolamer_pre_pos_).atom_index("N"), triazolamer_pre_pos_ );
		AtomID aidCA( pose.residue(triazolamer_pre_pos_).atom_index("CA"), triazolamer_pre_pos_ );
		AtomID aidC( pose.residue(triazolamer_pre_pos_).atom_index("C"), triazolamer_pre_pos_ );

		pose.conformation().set_torsion_angle( aidCYP, aidN, aidCA, aidC, radians( offset ) );
	}


	TR << "current triazolamer_pre ("<< triazolamer_pre_pos_<< ") PHI: " << pose.phi(triazolamer_pre_pos_) << std::endl;
	pose.set_phi(triazolamer_pre_pos_, phi_angle_);
	TR << "new triazolamer_pre ("<< triazolamer_pre_pos_<< ") PHI: " << pose.phi(triazolamer_pre_pos_) << std::endl;
	//pose.dump_pdb( "rosetta_out_triazolamer_pre.pdb" );
	TR << "current triazolamer_pre ("<< triazolamer_pre_pos_<< ") PSI: " << pose.psi(triazolamer_pre_pos_) << std::endl;
	pose.set_psi(triazolamer_pre_pos_, psi_angle_);
	TR << "new triazolamer_pre ("<< triazolamer_pre_pos_<< ") PSI: " << pose.psi(triazolamer_pre_pos_) << std::endl;
	//pose.dump_pdb( "rosetta_out_triazolamer_pre_pos_t.pdb" );

	update_hydrogens_( pose );

}

std::string
TriazolamerMover::get_name() const {
	return "TriazolamerMover";
}

/// @brief
TriazolamerMover::TriazolamerMover(
	core::Size triazolamer_seq_position
): Mover(), triazolamer_pre_pos_(triazolamer_seq_position), triazolamer_post_pos_(triazolamer_seq_position+1)
{
	Mover::type( "TriazolamerMover" );
}

TriazolamerMover::TriazolamerMover(
	core::Size triazolamer_seq_position,
	core::Real phi_angle,
	core::Real psi_angle
): Mover(), triazolamer_pre_pos_(triazolamer_seq_position), triazolamer_post_pos_(triazolamer_seq_position+1), phi_angle_(phi_angle), psi_angle_(psi_angle)
{
	Mover::type( "TriazolamerMover" );

}

TriazolamerMover::~TriazolamerMover(){}

void TriazolamerMover::update_hydrogens_( core::pose::Pose & pose )
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	AtomID aidVZP( pose.residue(triazolamer_pre_pos_).atom_index("VZP"), triazolamer_pre_pos_);
	AtomID aidCYP( pose.residue(triazolamer_pre_pos_).atom_index("CYP"), triazolamer_pre_pos_);
	AtomID aidN( pose.residue(triazolamer_pre_pos_).atom_index("N"), triazolamer_pre_pos_);
	AtomID aidCA( pose.residue(triazolamer_pre_pos_).atom_index("CA"), triazolamer_pre_pos_);
	AtomID aid1HYP( pose.residue(triazolamer_pre_pos_).atom_index("1HYP"), triazolamer_pre_pos_);
	AtomID aidVYP( pose.residue(triazolamer_post_pos_).atom_index("VYP"), triazolamer_post_pos_);
	AtomID aidCZP( pose.residue(triazolamer_post_pos_).atom_index("CZP"), triazolamer_post_pos_);
	AtomID aidN_Z( pose.residue(triazolamer_post_pos_).atom_index("N"), triazolamer_post_pos_);
	AtomID aidCA_Z( pose.residue(triazolamer_post_pos_).atom_index("CA"), triazolamer_post_pos_);
	AtomID aid1HZP( pose.residue(triazolamer_post_pos_).atom_index("1HZP"), triazolamer_post_pos_);

	//kdrew: definitions for xyz coordinates
	Vector const& CA_xyz ( pose.residue(triazolamer_pre_pos_).xyz("CA") );
	Vector const& N_xyz ( pose.residue(triazolamer_pre_pos_).xyz("N") );
	Vector const& CYP_xyz ( pose.residue(triazolamer_pre_pos_).xyz("CYP") );
	Vector const& CAZ_xyz ( pose.residue(triazolamer_post_pos_).xyz("CA") );
	Vector const& NZ_xyz ( pose.residue(triazolamer_post_pos_).xyz("N") );
	Vector const& CZP_xyz ( pose.residue(triazolamer_post_pos_).xyz("CZP") );


	//kdrew: VZP and VYP are relative to the 1Hs, so move the 1Hs and VZP|VYP and 2H moves with it
	//kdrew: calculate the difference between the virtual atom VZP torsion and the torsion involving the real CZP atom (on the post residue)
	Real VZP_torsion_correction = numeric::dihedral_degrees( CZP_xyz, CYP_xyz, N_xyz, CA_xyz ) - degrees( pose.conformation().torsion_angle( aidVZP, aidCYP, aidN, aidCA ) );
	//kdrew: change the 1HYP torsion by the corrected amount
	Real torsion_1HYP = degrees(pose.conformation().torsion_angle(aid1HYP,aidCYP,aidN,aidCA)) + VZP_torsion_correction;
	pose.conformation().set_torsion_angle(aid1HYP,aidCYP,aidN,aidCA,radians(torsion_1HYP));
	//pose.dump_pdb( "rosetta_out_triazolamer_pre_mvH.pdb" );

	//kdrew: calculate the difference between the virtual atom VYP torsion and the torsion involving the real CYP atom (on the pre residue)
	Real VYP_torsion_correction = numeric::dihedral_degrees( CYP_xyz, CZP_xyz, NZ_xyz, CAZ_xyz ) - degrees( pose.conformation().torsion_angle( aidVYP, aidCZP, aidN_Z, aidCA_Z ));
	//kdrew: change the 1HZP torsion by the corrected amount
	Real torsion_1HZP = degrees(pose.conformation().torsion_angle(aid1HZP,aidCZP,aidN_Z,aidCA_Z)) + VYP_torsion_correction;
	pose.conformation().set_torsion_angle(aid1HZP,aidCZP,aidN_Z,aidCA_Z,radians(torsion_1HZP));
	//pose.dump_pdb( "rosetta_out_triazolamer_pre_pos_t_mvH.pdb" );

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

}//triazolamer
}//simple_moves
}//protocols

