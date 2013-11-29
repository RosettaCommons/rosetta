// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/oop/OopPatcher.cc
/// @brief OopPatcher methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/oop/OopPatcher.hh>
// Package Headers

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/ncbb/util.hh>
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

static basic::Tracer TR( "protocols.simple_moves.oop.OopPatcher" );


using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace oop {

void OopPatcher::apply( core::pose::Pose & pose )
{
	TR<< "patching residues" <<std::endl;
	
	//kdrew: an oop pre position cannot be last position
	runtime_assert_msg ( oop_pre_pos_ != pose.total_residue(), "beginning of oop cannot be last residue" );
	//kdrew: an oop post position cannot be first position
	runtime_assert ( oop_post_pos_ != 1 );

	chemical::ResidueTypeSetCAP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	std::string const pre_base_name( core::chemical::residue_type_base_name( pose.residue_type( oop_pre_pos_ ) ) );
	std::string const post_base_name( core::chemical::residue_type_base_name( pose.residue_type( oop_post_pos_ ) ) );
	TR << "pre restype basename: " << pre_base_name << std::endl;
	TR << "post restype basename: " << post_base_name << std::endl;

	//kdrew: check for proline
	if ( pre_base_name == "PRO" || pre_base_name == "DPRO" ||
		 post_base_name == "PRO" || post_base_name == "DPRO" )
	{
    	utility_exit_with_message("Cannot patch proline");
	}
	if ( pose.residue(oop_pre_pos_).has_variant_type(chemical::OOP_POST) == 1) 
	{
    	utility_exit_with_message("Cannot patch OOP_PRE on an OOP_POST");
	}
	if ( pose.residue(oop_post_pos_).has_variant_type(chemical::OOP_PRE) == 1) 
	{
    	utility_exit_with_message("Cannot patch OOP_POST on an OOP_PRE");
	}

	//kdrew: check if already patched
	if ( pose.residue(oop_pre_pos_).has_variant_type(chemical::OOP_PRE) != 1)
	{
		TR<< "patching pre" <<std::endl;

		//kdrew: get base residue type
		chemical::ResidueType const & pre_base_type = pose.residue(oop_pre_pos_).type();
		TR<< pre_base_type.name() << std::endl;

		//kdrew: add variant
		conformation::Residue replace_res_pre( restype_set->get_residue_type_with_variant_added(pre_base_type, chemical::OOP_PRE), true );
		TR<< replace_res_pre.name() << std::endl;

		replace_res_pre.set_all_chi(pose.residue(oop_pre_pos_).chi());
		//replace_res_pre.mainchain_torsions(pose.residue(oop_pre_pos_).mainchain_torsions());

		pose.replace_residue( oop_pre_pos_, replace_res_pre, true );
		conformation::idealize_position( oop_pre_pos_, pose.conformation() );
		//pose.dump_pdb( "rosetta_out_oop_post_patch.pdb" );

	}// if pre
	if ( pose.residue(oop_post_pos_).has_variant_type(chemical::OOP_POST) != 1 ) {
		TR<< "patching post" <<std::endl;
		//kdrew: get base residue type
		chemical::ResidueType const & post_base_type = pose.residue(oop_post_pos_).type();
		TR<< post_base_type.name() << std::endl;

		//kdrew: add variant
		conformation::Residue replace_res_post( restype_set->get_residue_type_with_variant_added(post_base_type, chemical::OOP_POST), true );

		replace_res_post.set_all_chi(pose.residue(oop_post_pos_).chi());
		//replace_res_post.mainchain_torsions(pose.residue(oop_post_pos_).mainchain_torsions());

		pose.replace_residue( oop_post_pos_, replace_res_post, true );
		conformation::idealize_position( oop_post_pos_, pose.conformation() );
		//pose.dump_pdb( "rosetta_out_oop_pre_pos_t_patch.pdb" );

	}// if post

	core::pose::ncbb::add_oop_constraint( pose, oop_pre_pos_ );

	//kdrew: need to do all at once at the end because occasionally will get error:
	//kdrew: Unable to handle change in the number of residue connections in the presence of pseudobonds!
	//pose.conformation().detect_bonds();
	//pose.conformation().detect_pseudobonds();
	//pose.conformation().update_polymeric_connection( oop_pre_pos_ );
	//pose.conformation().update_polymeric_connection( oop_post_pos_ );

}


std::string
OopPatcher::get_name() const {
	return "OopPatcher";
}

///@brief
OopPatcher::OopPatcher(
		core::Size oop_seq_position
	): Mover(), oop_pre_pos_(oop_seq_position), oop_post_pos_(oop_seq_position+1)
{
	Mover::type( "OopPatcher" );

}

OopPatcher::~OopPatcher(){}

}//oop
}//simple_moves
}//protocols

