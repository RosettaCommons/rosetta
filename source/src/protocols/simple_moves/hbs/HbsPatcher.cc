// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/hbs/HbsPatcher.cc
/// @brief HbsPatcher methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
// Package Headers

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
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

static basic::Tracer TR( "protocols.simple_moves.hbs.HbsPatcher" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace hbs {

void add_hbs_constraint( core::pose::Pose & pose, core::Size hbs_pre_position )
{
	add_hbs_constraint( pose, hbs_pre_position, 1.52, 0.05 );
}
void add_hbs_constraint( core::pose::Pose & pose, core::Size hbs_pre_position, core::Real distance, core::Real std )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	//kdrew: add constraint
	core::scoring::func::HarmonicFuncOP harm_func  (new core::scoring::func::HarmonicFunc( distance, std ) );
	core::scoring::func::HarmonicFuncOP harm_func_0  (new core::scoring::func::HarmonicFunc( 0, std ) );
	core::scoring::func::CircularHarmonicFuncOP ang_func  (new core::scoring::func::CircularHarmonicFunc( 3.14159*2/3, 0.02 ) );
	core::scoring::func::CircularHarmonicFuncOP dih_func  (new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.02 ) );

	AtomID aidCYH( pose.residue( hbs_pre_position ).atom_index("CYH"), hbs_pre_position );
	AtomID aidCZH( pose.residue( hbs_pre_position+2 ).atom_index("CZH"), hbs_pre_position+2 );
	AtomID aidVZH( pose.residue( hbs_pre_position ).atom_index("VZH"), hbs_pre_position );
	AtomID aidVYH( pose.residue( hbs_pre_position+2 ).atom_index("VYH"), hbs_pre_position+2 );
	AtomID aidCY2( pose.residue( hbs_pre_position ).atom_index("CY2"), hbs_pre_position );
	AtomID aidCY1( pose.residue( hbs_pre_position ).atom_index("CY1"), hbs_pre_position );

	ConstraintCOP atompair = new AtomPairConstraint( aidCYH, aidCZH, harm_func );
	ConstraintCOP atompair2 = new AtomPairConstraint( aidCYH, aidVYH, harm_func_0 );
	ConstraintCOP atompair3 = new AtomPairConstraint( aidCZH, aidVZH, harm_func_0 );
	ConstraintCOP angle = new AngleConstraint( aidCZH, aidCYH, aidCY2, ang_func );
	ConstraintCOP dihedral = new DihedralConstraint( aidCZH, aidCYH, aidCY2, aidCY1, dih_func );

	pose.add_constraint( atompair );
	pose.add_constraint( atompair2 );
	pose.add_constraint( atompair3 );
	pose.add_constraint( angle );
	pose.add_constraint( dihedral );

	TR << "added atom pair constraint to hbs with distance: " << distance << " and std: "<< std << std::endl;
	TR << "and atom pair constraints with the virtual atoms" << std::endl;

}

void HbsPatcher::apply( core::pose::Pose & pose )
{
	TR<< "patching residues" <<std::endl;

	//kdrew: an hbs pre position cannot be last position
	runtime_assert_msg ( hbs_pre_pos_ != pose.total_residue(), "beginning of hbs cannot be last residue" );
	// I believe that this should be terminal, but since we're manually cutting off later residues
	// I bet we can't require it that way, so we have to trust hbs creator not to be dumb.
	//kdrew: an hbs post position cannot be first position
	runtime_assert ( hbs_post_pos_ != 1 );

	chemical::ResidueTypeSetCAP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	std::string const pre_base_name( core::chemical::residue_type_base_name( pose.residue_type( hbs_pre_pos_ ) ) );
	std::string const post_base_name( core::chemical::residue_type_base_name( pose.residue_type( hbs_post_pos_ ) ) );
	TR << "pre restype basename: " << pre_base_name << std::endl;
	TR << "post restype basename: " << post_base_name << std::endl;

	//kdrew: check for proline
	if ( pre_base_name == "PRO" || pre_base_name == "DPRO" ||
		 post_base_name == "PRO" || post_base_name == "DPRO" )
	{
    	utility_exit_with_message("Cannot patch proline");
	}
	if ( pose.residue(hbs_pre_pos_).has_variant_type(chemical::HBS_POST) == 1)
	{
    	utility_exit_with_message("Cannot patch HBS_PRE on an HBS_POST");
	}
	if ( pose.residue(hbs_post_pos_).has_variant_type(chemical::HBS_PRE) == 1)
	{
    	utility_exit_with_message("Cannot patch HBS_POST on an HBS_PRE");
	}

	//kdrew: check if already patched
	if ( pose.residue(hbs_pre_pos_).has_variant_type(chemical::HBS_PRE) != 1)
	{
		TR<< "patching pre" <<std::endl;

		//kdrew: get base residue type
		chemical::ResidueType const & pre_base_type = pose.residue(hbs_pre_pos_).type();
		TR<< pre_base_type.name() << std::endl;

		//kdrew: add variant
		conformation::Residue replace_res_pre( restype_set->get_residue_type_with_variant_added(pre_base_type, chemical::HBS_PRE), true );
		TR<< replace_res_pre.name() << std::endl;

		replace_res_pre.set_all_chi(pose.residue(hbs_pre_pos_).chi());
		//replace_res_pre.mainchain_torsions(pose.residue(hbs_pre_pos_).mainchain_torsions());

		pose.replace_residue( hbs_pre_pos_, replace_res_pre, true );
		conformation::idealize_position( hbs_pre_pos_, pose.conformation() );
		//pose.dump_pdb( "rosetta_out_hbs_post_patch.pdb" );

	}// if pre
	if ( pose.residue(hbs_post_pos_).has_variant_type(chemical::HBS_POST) != 1 ) {
		TR<< "patching post" <<std::endl;
		//kdrew: get base residue type
		chemical::ResidueType const & post_base_type = pose.residue(hbs_post_pos_).type();
		TR<< post_base_type.name() << std::endl;

		//kdrew: add variant
		conformation::Residue replace_res_post( restype_set->get_residue_type_with_variant_added(post_base_type, chemical::HBS_POST), true );

		replace_res_post.set_all_chi(pose.residue(hbs_post_pos_).chi());
		//replace_res_post.mainchain_torsions(pose.residue(hbs_post_pos_).mainchain_torsions());

		pose.replace_residue( hbs_post_pos_, replace_res_post, true );
		conformation::idealize_position( hbs_post_pos_, pose.conformation() );
		//pose.dump_pdb( "rosetta_out_hbs_pre_pos_t_patch.pdb" );

	}// if post

	add_hbs_constraint( pose, hbs_pre_pos_ );

	//kdrew: need to do all at once at the end because occasionally will get error:
	//kdrew: Unable to handle change in the number of residue connections in the presence of pseudobonds!
	//pose.conformation().detect_bonds();
	//pose.conformation().detect_pseudobonds();
	//pose.conformation().update_polymeric_connection( hbs_pre_pos_ );
	//pose.conformation().update_polymeric_connection( hbs_post_pos_ );

}


std::string
HbsPatcher::get_name() const {
	return "HbsPatcher";
}

///@brief
HbsPatcher::HbsPatcher(
		core::Size hbs_seq_position
	): Mover(), hbs_pre_pos_(hbs_seq_position), hbs_post_pos_(hbs_seq_position+2)
{
	Mover::type( "HbsPatcher" );

}

HbsPatcher::~HbsPatcher(){}

}//hbs
}//simple_moves
}//protocols
