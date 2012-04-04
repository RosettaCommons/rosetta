// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/oop/OopPuckMover.cc
/// @brief OopPuckMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/oop/OopPuckMover.fwd.hh>
#include <protocols/simple_moves/oop/OopPuckMover.hh>
#include <protocols/simple_moves/oop/OopMover.hh>
#include <protocols/simple_moves/oop/OopPatcher.hh>
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

static numeric::random::RandomGenerator RG(956732);
static basic::Tracer TR( "protocols.simple_moves.oop.OopPuckMover" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

//kdrew: defining constants
static const Real OOP_PUCK_PLUS_PHI = -131.27;
static const Real OOP_PUCK_PLUS_PSI = -10.38;

static const Real OOP_PUCK_MINUS_PHI = -147.05;
static const Real OOP_PUCK_MINUS_PSI = -36.90;


/*
void oop_patch ( core::pose::Pose & pose, core::Size oop_pre_pos )
{

	core::Size oop_post_pos = oop_pre_pos + 1;

	TR<< "patching residues" <<std::endl;
	chemical::ResidueTypeSetCAP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	//kdrew: check if already patched
	if ( pose.residue(oop_pre_pos).has_variant_type(chemical::OOP_PRE) != 1) 
	{
		TR<< "patching pre" <<std::endl;

		//kdrew: get base residue type
		chemical::ResidueType pre_base_type = pose.residue(oop_pre_pos).type();
		TR<< pre_base_type.name() << std::endl;

		//kdrew: add variant
		conformation::Residue replace_res_pre( restype_set->get_residue_type_with_variant_added(pre_base_type, chemical::OOP_PRE), true );

		replace_res_pre.set_all_chi(pose.residue(oop_pre_pos).chi());
		//replace_res_pre.mainchain_torsions(pose.residue(oop_pre_pos).mainchain_torsions());
		

		pose.replace_residue( oop_pre_pos, replace_res_pre, true );
		conformation::idealize_position( oop_pre_pos, pose.conformation() );
		//pose.dump_pdb( "rosetta_out_oop_post_patch.pdb" );

	}// if pre
	if ( pose.residue(oop_post_pos).has_variant_type(chemical::OOP_POST) != 1 ) {
		TR<< "patching post" <<std::endl;
		//kdrew: get base residue type
		chemical::ResidueType post_base_type = pose.residue(oop_post_pos).type();
		TR<< post_base_type.name() << std::endl;

		//kdrew: add variant
		conformation::Residue replace_res_post( restype_set->get_residue_type_with_variant_added(post_base_type, chemical::OOP_POST), true );

		replace_res_post.set_all_chi(pose.residue(oop_post_pos).chi());
		//replace_res_post.mainchain_torsions(pose.residue(oop_post_pos).mainchain_torsions());

		pose.replace_residue( oop_post_pos, replace_res_post, true );
		conformation::idealize_position( oop_post_pos, pose.conformation() );
		//pose.dump_pdb( "rosetta_out_oop_pre_post_patch.pdb" );

	}// if post

}
*/

/*
kdrew: the helper function moves a single oop with the given move
pose: pose to make change to
oop_pre_pos: position of first residue of oop to make change to
phi_angle: value to change phi angle to
psi_angle: value to change psi angle to
void oop_puck_mover_helper( core::pose::Pose & pose, core::Size oop_pre_pos, Real phi_angle, Real psi_angle)
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	core::Size oop_post_pos = oop_pre_pos + 1;

	pose.set_phi(oop_pre_pos, phi_angle);
	//pose.dump_pdb( "rosetta_out_oop_pre.pdb" );
	pose.set_psi(oop_pre_pos, psi_angle);
	//pose.dump_pdb( "rosetta_out_oop_pre_post.pdb" );

	AtomID aidVZP( pose.residue(oop_pre_pos).atom_index("VZP"), oop_pre_pos);
	AtomID aidCYP( pose.residue(oop_pre_pos).atom_index("CYP"), oop_pre_pos);
	AtomID aidN( pose.residue(oop_pre_pos).atom_index("N"), oop_pre_pos);
	AtomID aidCA( pose.residue(oop_pre_pos).atom_index("CA"), oop_pre_pos);
	AtomID aid1HYP( pose.residue(oop_pre_pos).atom_index("1HYP"), oop_pre_pos);
	AtomID aidVYP( pose.residue(oop_post_pos).atom_index("VYP"), oop_post_pos);
	AtomID aidCZP( pose.residue(oop_post_pos).atom_index("CZP"), oop_post_pos);
	AtomID aidN_Z( pose.residue(oop_post_pos).atom_index("N"), oop_post_pos);
	AtomID aidCA_Z( pose.residue(oop_post_pos).atom_index("CA"), oop_post_pos);
	AtomID aid1HZP( pose.residue(oop_post_pos).atom_index("1HZP"), oop_post_pos);

	//kdrew: definitions for xyz coordinates 
	Vector const& CA_xyz ( pose.residue(oop_pre_pos).xyz("CA") );
	Vector const& N_xyz ( pose.residue(oop_pre_pos).xyz("N") );
	Vector const& CYP_xyz ( pose.residue(oop_pre_pos).xyz("CYP") );
	Vector const& CAZ_xyz ( pose.residue(oop_post_pos).xyz("CA") );
	Vector const& NZ_xyz ( pose.residue(oop_post_pos).xyz("N") );
	Vector const& CZP_xyz ( pose.residue(oop_post_pos).xyz("CZP") );


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
	//pose.dump_pdb( "rosetta_out_oop_pre_post_mvH.pdb" );

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

	//pose.dump_pdb( "rosetta_moved_puck_helper.pdb"  );

}
*/

namespace protocols {
namespace simple_moves {
namespace oop {
/*
///@details
void OopPuckMover::apply( core::pose::Pose & pose ){
	
	TR<< "in OopPuckMover::apply" << std::endl;
	//kdrew: for all positions in oop_seq_positions_, input assertion check
	for(Size i = 1; i <= oop_seq_positions_.size(); i++)
	{
		Size oop_pre_pos = oop_seq_positions_[i];
		Size oop_post_pos = oop_pre_pos+1;                   
		TR<< "oop_pre_pos:" << oop_pre_pos << " oop_post_pos:" << oop_post_pos << std::endl;

		//kdrew: patch will add variant types later
		if ( !patch_ )
		{
			runtime_assert ( pose.residue(oop_pre_pos).has_variant_type(chemical::OOP_PRE) == 1) ;
			runtime_assert ( pose.residue(oop_post_pos).has_variant_type(chemical::OOP_POST) == 1) ;
		}
		//kdrew: an oop pre position cannot be last position
		runtime_assert ( oop_pre_pos != pose.total_residue() );
		//kdrew: an oop post position cannot be first position
		runtime_assert ( oop_post_pos != 1 );

	}//for


	//kdrew: randomly pick positions and puckers to move
	if( random_ )
	{
		//kdrew: randomly choose position from oop_seq_positions
		core::Size random_pos = oop_seq_positions_[int(RG.uniform()*oop_seq_positions_.size())+1];

		//kdrew: randomly choose conformation up, down or small angle move
		std::string random_pucker = available_moves_[int(RG.uniform()*available_moves_.size())+1];

		//kdrew: add oop patch if residue does not have it
		if ( patch_ ) {
			oop::OopPatcherOP oop_patcher (new oop::OopPatcher( random_pos ) );
			oop_patcher->apply( pose );
		}

		runtime_assert ( random_pucker == "OOP_PUCK_PLUS" || random_pucker == "OOP_PUCK_MINUS" || random_pucker == "OOP_PUCK_SMALL" );
		TR << random_pucker <<std::endl;

        oop::OopMoverOP oop_mover ( new oop::OopMover( random_pos ) );
        oop::OopPuckPlusMoverOP oop_plus_mover ( new oop::OopPuckPlusMover( random_pos ) );
        oop::OopPuckMinusMoverOP oop_minus_mover ( new oop::OopPuckMinusMover( random_pos ) );

		if ( random_pucker == "OOP_PUCK_PLUS" ) 
		{
			oop_plus_mover->apply( pose );
			//oop_puck_mover_helper(pose, random_pos, OOP_PUCK_PLUS_PHI, OOP_PUCK_PLUS_PSI);
		}
		else if (random_pucker == "OOP_PUCK_MINUS" )
		{
			oop_minus_mover->apply( pose );
			//oop_puck_mover_helper(pose, random_pos, OOP_PUCK_MINUS_PHI, OOP_PUCK_MINUS_PSI);
		}

		else if (random_pucker == "OOP_PUCK_SMALL" )
		{
			Real small_angle = max_small_angle_/2.0; ///< this is max_angle/2, which is the deviation from the angle input
			Real phi_angle = basic::periodic_range( pose.phi( random_pos ) - small_angle + RG.uniform() * max_small_angle_, 360.0 );
			Real psi_angle = basic::periodic_range( pose.psi( random_pos ) - small_angle + RG.uniform() * max_small_angle_, 360.0 );
			oop_mover->set_phi( phi_angle );
			oop_mover->set_psi( psi_angle );
			oop_mover->apply( pose );
			//oop_puck_mover_helper( pose, random_pos, phi_angle, psi_angle );
		}

	}//if
	else //change all positions
	{
		for(Size i = 1; i <= oop_seq_positions_.size(); i++)
		{
			if( available_moves_.size() > 0 )
			{
				//kdrew: add oop patch if residue does not have it
				if ( patch_ ) {
					oop::OopPatcherOP oop_patcher (new oop::OopPatcher( oop_seq_positions_[i] ) );
					oop_patcher->apply( pose );
				}

				oop::OopMoverOP oop_mover ( new oop::OopMover( oop_seq_positions_[i] ) );

				//kdrew: this is ugly code duplication from above
				runtime_assert ( available_moves_[1] == "OOP_PUCK_PLUS" || available_moves_[1] == "OOP_PUCK_MINUS" || available_moves_[1] == "OOP_PUCK_SMALL" );
				TR << available_moves_[1] <<std::endl;

				if ( available_moves_[1] == "OOP_PUCK_PLUS" ) 
				{
					oop_mover->set_phi( OOP_PUCK_PLUS_PHI );
					oop_mover->set_psi( OOP_PUCK_PLUS_PSI );
					oop_mover->apply( pose );
					//oop_puck_mover_helper(pose, oop_seq_positions_[i], OOP_PUCK_PLUS_PHI, OOP_PUCK_PLUS_PSI);
				}
				else if (available_moves_[1] == "OOP_PUCK_MINUS" )
				{
					oop_mover->set_phi( OOP_PUCK_MINUS_PHI );
					oop_mover->set_psi( OOP_PUCK_MINUS_PSI );
					oop_mover->apply( pose );
					//oop_puck_mover_helper(pose, oop_seq_positions_[i], OOP_PUCK_MINUS_PHI, OOP_PUCK_MINUS_PSI);
				}

				else if (available_moves_[1] == "OOP_PUCK_SMALL" )
				{
					Real small_angle = max_small_angle_/2.0; ///< this is max_angle/2, which is the deviation from the angle input
					Real phi_angle = basic::periodic_range( pose.phi( oop_seq_positions_[i] ) - small_angle + RG.uniform() * max_small_angle_, 360.0 );
					Real psi_angle = basic::periodic_range( pose.psi( oop_seq_positions_[i] ) - small_angle + RG.uniform() * max_small_angle_, 360.0 );
					oop_mover->set_phi( phi_angle );
					oop_mover->set_psi( psi_angle );
					oop_mover->apply( pose );
					//oop_puck_mover_helper( pose, oop_seq_positions_[i], phi_angle, psi_angle );
				}
			}

		}//for

	}//else

}//apply

std::string
OopPuckMover::get_name() const {
	return "OopPuckMover";
}

///@brief
OopPuckMover::OopPuckMover(
) : Mover()
{
	Mover::type( "OopPuckMover" );
}

OopPuckMover::OopPuckMover( 
		utility::vector1< core::Size > oop_seq_positions, 
		bool oop_puck_plus, 
		bool oop_puck_minus,
		bool random,
		bool patch
	): Mover(), oop_seq_positions_(oop_seq_positions), oop_puck_plus_(oop_puck_plus), oop_puck_minus_(oop_puck_minus), random_(random), patch_(patch), oop_puck_small_(false)
{
	Mover::type( "OopPuckMover" );

	if (oop_puck_plus_ )
	{
    	available_moves_.push_back("OOP_PUCK_PLUS");	
	}
	if (oop_puck_minus_ )
	{
    	available_moves_.push_back("OOP_PUCK_MINUS");	
	}
}

OopPuckMover::OopPuckMover( 
		utility::vector1< core::Size > oop_seq_positions, 
		bool oop_puck_plus, 
		bool oop_puck_minus,
		bool random,
		bool patch,
		bool oop_puck_small,
		Real max_small_angle
	): Mover(), oop_seq_positions_(oop_seq_positions), oop_puck_plus_(oop_puck_plus), oop_puck_minus_(oop_puck_minus), random_(random), patch_(patch), oop_puck_small_(oop_puck_small), max_small_angle_(max_small_angle)
{
	Mover::type( "OopPuckMover" );

	if (oop_puck_plus_ )
	{
    	available_moves_.push_back("OOP_PUCK_PLUS");	
	}
	if (oop_puck_minus_ )
	{
    	available_moves_.push_back("OOP_PUCK_MINUS");	
	}
	if ( oop_puck_small_ )
	{
    	available_moves_.push_back("OOP_PUCK_SMALL");	
	}
}

OopPuckMover::OopPuckMover( utility::vector1< core::Size > oop_seq_positions ): Mover(), oop_seq_positions_(oop_seq_positions), oop_puck_plus_(true), oop_puck_minus_(true), random_(true), patch_(false), oop_puck_small_(false), max_small_angle_(0.0)
{
	Mover::type( "OopPuckMover" );
	available_moves_.push_back("OOP_PUCK_PLUS");	
	available_moves_.push_back("OOP_PUCK_MINUS");	
}

OopPuckMover::~OopPuckMover(){}
*/


OopPuckPlusMover::OopPuckPlusMover( core::Size oop_seq_position ): OopMover( oop_seq_position, OOP_PUCK_PLUS_PHI, OOP_PUCK_PLUS_PSI )
{
	OopMover::type( "OopPuckPlusMover" );
}

OopPuckPlusMover::~OopPuckPlusMover(){}

std::string
OopPuckPlusMover::get_name() const {
	return "OopPuckPlusMover";
}

OopPuckMinusMover::OopPuckMinusMover( core::Size oop_seq_position ): OopMover( oop_seq_position, OOP_PUCK_MINUS_PHI, OOP_PUCK_MINUS_PSI )
{
	OopMover::type( "OopPuckMinusMover" );
}

OopPuckMinusMover::~OopPuckMinusMover(){}

std::string
OopPuckMinusMover::get_name() const {
	return "OopPuckMinusMover";
}


}//oop
}//simple_moves
}//protocols

