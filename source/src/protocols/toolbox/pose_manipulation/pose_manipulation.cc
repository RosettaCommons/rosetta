// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/pose_manipulation.cc
/// @brief some general functions to manipulate poses. feel free to add your own
/// @brief if you add your own, please mention your name:)
/// @author Florian Richter, floric@u.washington.edu


// Unit headers
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/util.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/id/AtomID_Map.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh> //typeset swapping

#include <protocols/moves/MonteCarlo.hh>

#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/Loops.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <core/chemical/VariantType.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>


namespace protocols {
namespace toolbox {
namespace pose_manipulation{

static thread_local basic::Tracer TR( "protocols.toolbox.pose_manipulation" );
static thread_local basic::Tracer TR_DI( "protocols.toolbox.pose_manipulation.insert_pose_into_pose" );
using basic::T;
using basic::Error;
using basic::Warning;
using core::chemical::ResidueType;

void
construct_poly_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
)
{
	construct_poly_uniq_restype_pose( pose, positions, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA"), keep_pro, keep_gly, keep_disulfide_cys);
}


void
construct_poly_uniq_restype_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	ResidueType const & restype,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
)
{
	using namespace core;

	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	conformation::Residue const replace_res( restype, true );

	for( utility::vector1< Size >::const_iterator pos_it = positions.begin();
			 pos_it != positions.end(); ++pos_it ) {

		chemical::ResidueType const * cur_restype = & pose.residue_type( *pos_it );


		if( ( keep_pro && ( cur_restype->aa() == chemical::aa_pro ) )
				||( keep_gly && ( cur_restype->aa() == chemical::aa_gly ) )
				||( keep_disulfide_cys && ( cur_restype->aa() == chemical::aa_cys ) && cur_restype->has_variant_type( chemical::DISULFIDE ) ) )
		{
			continue;
		}

		utility::vector1< std::string > current_variants;

		if ( TR.Debug.visible() ) {
			TR.Debug << "replacing: " << *pos_it << std::endl;
		}

		// either we don't want to keep disulfide cys or the current restype is not cys
		// so ignore disulfide variant type
		if ( ! variants_match_with_exceptions( *cur_restype, replace_res.type(), utility::vector1< core::chemical::VariantType >( 1, core::chemical::DISULFIDE ) ) ) {
				current_variants = cur_restype->properties().get_list_of_variants();
				chemical::ResidueTypeCOP var_replace_type = replace_res.type().get_self_ptr();

			for ( core::Size var = 1; var <= current_variants.size(); ++var ) {
				if ( current_variants[ var ] != "DISULFIDE" ) {
					var_replace_type = restype_set->get_residue_type_with_variant_added( * var_replace_type,
							core::chemical::ResidueProperties::get_variant_from_string( current_variants[ var ] ) ).get_self_ptr();
				}
			}

			//runtime_assert( var_replace_type->name3() == "ALA" );
			conformation::Residue const var_replace_res( *var_replace_type, true );

			pose.replace_residue( *pos_it, var_replace_res, true );
		} else {
			pose.replace_residue( *pos_it, replace_res, true);
		}

	} //iterator over positions to replace

} // construct_poly_ala_pose function

void
construct_poly_XXX_pose(
	std::string const & aa,
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
)
{
	using namespace core;

	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	construct_poly_XXX_pose( aa, pose, positions, restype_set, keep_pro, keep_gly, keep_disulfide_cys );
}

void
construct_poly_XXX_pose(
	std::string const & aa,
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	core::chemical::ResidueTypeSetCOP restype_set,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
)
{
	using namespace core;
	
	conformation::Residue const replace_res( restype_set->name_map( aa ), true );

	for ( utility::vector1< Size >::const_iterator pos_it = positions.begin();
			 pos_it != positions.end(); ++pos_it ) {

		chemical::ResidueType const * cur_restype = & pose.residue_type( *pos_it );

		if( ( keep_pro && ( cur_restype->aa() == chemical::aa_pro ) )
				||( keep_gly && ( cur_restype->aa() == chemical::aa_gly ) )
				||( keep_disulfide_cys && ( cur_restype->aa() == chemical::aa_cys ) && cur_restype->has_variant_type( chemical::DISULFIDE ) ) )
		{
			continue;
		}
		utility::vector1< std::string > current_variants;

		if ( ! variants_match( *cur_restype, replace_res.type() ) ) {
			current_variants = cur_restype->properties().get_list_of_variants();
			chemical::ResidueTypeCOP var_replace_type = replace_res.type().get_self_ptr();

			for ( core::Size var = 1; var <= current_variants.size(); ++var ) {
				if ( ( cur_restype->has_variant_type( chemical::DISULFIDE ) && keep_disulfide_cys) ||
					( ! cur_restype->has_variant_type( chemical::DISULFIDE ) ) ) {
					var_replace_type = restype_set->get_residue_type_with_variant_added( * var_replace_type,
							chemical::ResidueProperties::get_variant_from_string( current_variants[ var ] ) ).get_self_ptr();
				}
			}

			runtime_assert( var_replace_type->name3() == aa );
			conformation::Residue const var_replace_res( *var_replace_type, true );
			pose.replace_residue( *pos_it, var_replace_res, true );
		} else {
			pose.replace_residue( *pos_it, replace_res, true);
		}
	} //iterator over positions to replace
} // construct_poly_XXX_pose function


void
remove_non_protein_residues(
	core::pose::Pose & pose
)
{

	bool residues_deleted(false);

	for( core::Size i = pose.total_residue(); i > 0 ; --i){

		if( ! pose.residue_type( i ).is_protein() ){
			pose.conformation().delete_residue_slow( i );
			residues_deleted = true;
		}
	}

	if( residues_deleted ) pose.energies().clear();

} //remove_non_protein_residues


void
add_chainbreaks_according_to_jumps( core::pose::Pose & pose )
{
	for( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ){

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );

		if ( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) ) continue;

		if ( pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) continue;

		if ( !pose.residue_type( this_cutpoint ).is_protein() ) continue;
		if ( !pose.residue_type( this_cutpoint +1 ).is_protein() ) continue;

		if( !pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ){
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if( !pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ){
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this_cutpoint +1 );
		}
	}
}

void
add_chainbreaks_according_to_jumps( core::pose::Pose & pose, utility::vector1< core::Size > const& no_cutpoint_residues  )
{
	for( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ){

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );
		//exclude residue numbers in array
		if ( find( no_cutpoint_residues.begin(), no_cutpoint_residues.end(), this_cutpoint ) != no_cutpoint_residues.end() ) {
			continue;
		}

		if ( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) ) continue;

		if ( pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) continue;

		if( !pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ){
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if( !pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ){
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this_cutpoint +1 );
		}
	}
}


void
remove_chainbreaks_according_to_jumps( core::pose::Pose & pose )
{
	for( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ){

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );

		if( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ){
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if( pose.residue_type( this_cutpoint + 1).has_variant_type( core::chemical::CUTPOINT_UPPER ) ){
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this_cutpoint+1 );
		}
	}
}

core::Real
superimpose_pose_on_subset_CA(
	core::pose::Pose & pose,
	core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & positions,
	int const offset
)
{
	using namespace core::id;

	AtomID_Map< AtomID > atom_map;

	core::pose::initialize_atomid_map( atom_map, pose, BOGUS_ATOM_ID );

	for( utility::vector1< core::Size >::const_iterator res_it = positions.begin(); res_it != positions.end(); ++res_it){

		AtomID id1( pose.residue( *res_it + offset).atom_index("CA"), *res_it + offset );
		AtomID id2( ref_pose.residue( *res_it ).atom_index("CA"), *res_it );

		atom_map.set( id1, id2);

	} //iterator over residues

	return core::scoring::superimpose_pose( pose, ref_pose, atom_map );
} //superimpose_pose_on_subset_CA


} // namespace pose_manipulation
} //toolbox
} //protocols
