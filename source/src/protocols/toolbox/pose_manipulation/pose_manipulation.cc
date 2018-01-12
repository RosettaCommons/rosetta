// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/pose/util.hh>

#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/select/util.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <core/id/AtomID_Map.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh> //typeset swapping

#include <protocols/moves/MonteCarlo.hh>

#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/Loops.hh>

// Utility Headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <core/chemical/VariantType.hh>
#include <core/pose/variant_util.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>


namespace protocols {
namespace toolbox {
namespace pose_manipulation {

static basic::Tracer TR( "protocols.toolbox.pose_manipulation" );
static basic::Tracer TR_DI( "protocols.toolbox.pose_manipulation.insert_pose_into_pose" );
using basic::Error;
using basic::Warning;
using core::chemical::ResidueType;

void
construct_poly_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
) {
	construct_poly_uniq_restype_pose( pose, positions, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA"), keep_pro, keep_gly, keep_disulfide_cys);
}

/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @brief puts in D-ala residues at the positions specified in the 'positions' input array
void
construct_poly_d_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
) {
	construct_poly_uniq_restype_pose( pose, positions, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("DALA"), keep_pro, keep_gly, keep_disulfide_cys);
}

/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @brief puts in beta-3-ala residues at the positions specified in the 'positions' input array.
void
construct_poly_beta_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
) {
	construct_poly_uniq_restype_pose( pose, positions, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("B3A"), keep_pro, keep_gly, keep_disulfide_cys);
}

/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @brief puts in D-beta-3-ala residues at the positions specified in the 'positions' input array
void
construct_poly_d_beta_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
) {
	construct_poly_uniq_restype_pose( pose, positions, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("DB3A"), keep_pro, keep_gly, keep_disulfide_cys);
}

void
construct_poly_uniq_restype_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	ResidueType const & restype,
	bool const keep_pro,
	bool const keep_gly,
	bool const keep_disulfide_cys
)
{
	using namespace core;

	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	conformation::Residue const replace_res( restype, true );

	for ( unsigned long position : positions ) {

		if ( replace_res.type().is_alpha_aa() && !pose.residue_type( position ).is_alpha_aa() ) continue;
		if ( replace_res.type().is_beta_aa() && !pose.residue_type( position ).is_beta_aa() ) continue;

		chemical::ResidueType const & cur_restype = pose.residue_type( position );


		if ( ( keep_pro && ( cur_restype.aa() == chemical::aa_pro || cur_restype.aa() == chemical::aa_dpr || cur_restype.aa() == chemical::aa_b3p ) )
				||( keep_gly && ( cur_restype.aa() == chemical::aa_gly || cur_restype.aa() == chemical::aa_b3g ) )
				||( keep_disulfide_cys && cur_restype.has_variant_type( chemical::DISULFIDE ) ) ) {
			continue;
		}

		utility::vector1< std::string > current_variants;

		if ( TR.Debug.visible() ) {
			TR.Debug << "replacing: " << position << std::endl;
		}

		// either we don't want to keep disulfide cys or the current restype is not cys
		// so ignore disulfide variant type
		if ( ! variants_match_with_exceptions( cur_restype, replace_res.type(), utility::vector1< core::chemical::VariantType >( 1, core::chemical::DISULFIDE ) ) ) {
			current_variants = cur_restype.properties().get_list_of_variants();
			chemical::ResidueTypeCOP var_replace_type = replace_res.type_ptr();

			for ( core::Size var = 1; var <= current_variants.size(); ++var ) {
				if ( current_variants[ var ] != "DISULFIDE" ) {
					var_replace_type = restype_set->get_residue_type_with_variant_added( * var_replace_type,
						core::chemical::ResidueProperties::get_variant_from_string( current_variants[ var ] ) ).get_self_ptr();
				}
			}

			conformation::Residue const var_replace_res( *var_replace_type, true );

			pose.replace_residue( position, var_replace_res, true );
		} else {
			pose.replace_residue( position, replace_res, true);
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

	for ( unsigned long position : positions ) {

		chemical::ResidueType const & cur_restype = pose.residue_type( position );

		if ( ( keep_pro && ( cur_restype.aa() == chemical::aa_pro ) )
				||( keep_gly && ( cur_restype.aa() == chemical::aa_gly ) )
				||( keep_disulfide_cys && ( cur_restype.aa() == chemical::aa_cys ) && cur_restype.has_variant_type( chemical::DISULFIDE ) ) ) {
			continue;
		}
		utility::vector1< std::string > current_variants;

		if ( ! variants_match( cur_restype, replace_res.type() ) ) {
			current_variants = cur_restype.properties().get_list_of_variants();
			chemical::ResidueTypeCOP var_replace_type = replace_res.type_ptr();

			for ( core::Size var = 1; var <= current_variants.size(); ++var ) {
				if ( ( cur_restype.has_variant_type( chemical::DISULFIDE ) && keep_disulfide_cys) ||
						( ! cur_restype.has_variant_type( chemical::DISULFIDE ) ) ) {
					var_replace_type = restype_set->get_residue_type_with_variant_added( * var_replace_type,
						chemical::ResidueProperties::get_variant_from_string( current_variants[ var ] ) ).get_self_ptr();
				}
			}

			runtime_assert( var_replace_type->name3() == aa );
			conformation::Residue const var_replace_res( *var_replace_type, true );
			pose.replace_residue( position, var_replace_res, true );
		} else {
			pose.replace_residue( position, replace_res, true);
		}
	} //iterator over positions to replace
} // construct_poly_XXX_pose function


void
remove_non_protein_residues(
	core::pose::Pose & pose
)
{

	bool residues_deleted(false);

	for ( core::Size i = pose.size(); i > 0 ; --i ) {

		if ( ! pose.residue_type( i ).is_protein() ) {
			pose.conformation().delete_residue_slow( i );
			residues_deleted = true;
		}
	}

	if ( residues_deleted ) pose.energies().clear();

} //remove_non_protein_residues


void
add_chainbreaks_according_to_jumps( core::pose::Pose & pose )
{
	for ( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ) {

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );

		if ( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) ) continue;

		if ( pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) continue;

		if ( !pose.residue_type( this_cutpoint ).is_protein() ) continue;
		if ( !pose.residue_type( this_cutpoint +1 ).is_protein() ) continue;

		if ( !pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if ( !pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this_cutpoint +1 );
		}
	}
}

void
add_chainbreaks_according_to_jumps( core::pose::Pose & pose, utility::vector1< core::Size > const& no_cutpoint_residues  )
{
	for ( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ) {

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );
		//exclude residue numbers in array
		if ( find( no_cutpoint_residues.begin(), no_cutpoint_residues.end(), this_cutpoint ) != no_cutpoint_residues.end() ) {
			continue;
		}

		if ( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) ) continue;

		if ( pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) continue;

		if ( !pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if ( !pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this_cutpoint +1 );
		}
	}
}


void
remove_chainbreaks_according_to_jumps( core::pose::Pose & pose )
{
	for ( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ) {

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );

		if ( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if ( pose.residue_type( this_cutpoint + 1).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
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

	core::pose::initialize_atomid_map( atom_map, pose, AtomID::BOGUS_ATOM_ID() );

	for ( unsigned long position : positions ) {

		AtomID id1( pose.residue( position + offset).atom_index("CA"), position + offset );
		AtomID id2( ref_pose.residue( position ).atom_index("CA"), position );

		atom_map.set( id1, id2);

	} //iterator over residues

	return core::scoring::superimpose_pose( pose, ref_pose, atom_map );
} //superimpose_pose_on_subset_CA


void
repack_this_residue(
	core::Size seq_pos,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP scorefxn,
	bool include_current /* = true */,
	std::string name1s_if_design /* = "" */ ) {

	using namespace core::select::residue_selector;
	using namespace core::pack::task::operation;

	core::pack::task::TaskFactoryOP local_tf( new core::pack::task::TaskFactory() );

	local_tf->push_back( TaskOperationCOP( new InitializeFromCommandline() ) );
	if ( include_current ) {
		local_tf->push_back( TaskOperationCOP( new IncludeCurrent() ) );
	}
	if ( name1s_if_design.length() > 0 ) {
		RestrictAbsentCanonicalAASOP restrict_absent( new RestrictAbsentCanonicalAAS() );
		restrict_absent->keep_aas( name1s_if_design );
		restrict_absent->include_residue( 0 ); // 0 means apply to all
		local_tf->push_back( restrict_absent );
	} else {
		local_tf->push_back( TaskOperationCOP( new RestrictToRepacking() ) );
	}

	ResidueSelectorOP the_residue( new ResidueIndexSelector( utility::to_string( seq_pos ) ) );


	ResLvlTaskOperationOP prevent_repacking( new PreventRepackingRLT() );
	OperateOnResidueSubsetOP inv_subset( new OperateOnResidueSubset( prevent_repacking, the_residue, true) );
	local_tf->push_back( inv_subset );

	protocols::minimization_packing::PackRotamersMoverOP repack(
		new protocols::minimization_packing::PackRotamersMover( scorefxn ) );
	repack->task_factory( local_tf );
	repack->apply( pose );
}


void
repack_these_residues(
	core::select::residue_selector::ResidueSubset const & subset,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP scorefxn,
	bool include_current /* = true */,
	std::string name1s_if_design /* = "" */) {

	using namespace core::select::residue_selector;
	using namespace core::pack::task::operation;

	core::pack::task::TaskFactoryOP local_tf( new core::pack::task::TaskFactory() );

	local_tf->push_back( TaskOperationCOP( new InitializeFromCommandline() ) );
	if ( include_current ) {
		local_tf->push_back( TaskOperationCOP( new IncludeCurrent() ) );
	}
	if ( name1s_if_design.length() > 0 ) {
		RestrictAbsentCanonicalAASOP restrict_absent( new RestrictAbsentCanonicalAAS() );
		restrict_absent->keep_aas( name1s_if_design );
		restrict_absent->include_residue( 0 ); // 0 means apply to all
		local_tf->push_back( restrict_absent );
	} else {
		local_tf->push_back( TaskOperationCOP( new RestrictToRepacking() ) );
	}

	ResidueSelectorCOP these_residues( core::select::get_residue_selector_from_subset( subset ) );

	ResLvlTaskOperationOP prevent_repacking( new PreventRepackingRLT() );
	OperateOnResidueSubsetOP inv_subset( new OperateOnResidueSubset( prevent_repacking, these_residues, true) );
	local_tf->push_back( inv_subset );

	protocols::minimization_packing::PackRotamersMoverOP repack(
		new protocols::minimization_packing::PackRotamersMover( scorefxn ) );
	repack->task_factory( local_tf );
	repack->apply( pose );
}

void
rigid_body_move(
	numeric::xyzVector<core::Real> const & rotation_unit_vector,
	core::Real angle_deg,
	numeric::xyzVector<core::Real> const & translation_vector,
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	numeric::xyzVector<core::Real> center_of_rotation
	/* = numeric::xyzVector<core::Real>(std::numeric_limits<double>::quiet_NaN(), 0, 0) */) {

	numeric::xyzMatrix<core::Real> rotation = numeric::rotation_matrix( rotation_unit_vector,
		angle_deg * numeric::constants::r::pi / 180. );

	rigid_body_move( rotation, translation_vector, pose, subset, center_of_rotation );
}

void
rigid_body_move(
	numeric::xyzVector<core::Real> const & rotation_unit_vector,
	core::Real angle_deg,
	numeric::xyzVector<core::Real> const & translation_vector,
	core::Real translation_scalar,
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	numeric::xyzVector<core::Real> center_of_rotation
	/* = numeric::xyzVector<core::Real>(std::numeric_limits<double>::quiet_NaN(), 0, 0) */) {

	numeric::xyzMatrix<core::Real> rotation = numeric::rotation_matrix( rotation_unit_vector,
		angle_deg * numeric::constants::r::pi / 180. );

	numeric::xyzVector<core::Real> translate = translation_vector;
	translate *= translation_scalar;

	rigid_body_move( rotation, translate, pose, subset, center_of_rotation );
}


void
rigid_body_move(
	numeric::xyzMatrix<core::Real> rotation,
	numeric::xyzVector<core::Real> const & translation_vector,
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	numeric::xyzVector<core::Real> center_of_rotation
	/* = numeric::xyzVector<core::Real>(std::numeric_limits<double>::quiet_NaN(), 0, 0) */) {

	if ( std::isnan( center_of_rotation.x() ) || std::isnan( center_of_rotation.y() )
			|| std::isnan( center_of_rotation.z() ) ) {
		center_of_rotation = core::pose::center_of_mass( pose, subset );
	}

	utility::vector1<core::Size> seq_poss = core::select::get_residues_from_subset( subset );

	for ( core::Size seq_pos : seq_poss ) {
		core::conformation::Residue const & res( pose.residue( seq_pos ) );
		core::Size natoms = res.natoms();
		for ( core::Size j = 1; j <= natoms; ++j ) {
			numeric::xyzVector<core::Real> working = res.atom(j).xyz();
			working -= center_of_rotation;
			working = rotation*working;
			working += center_of_rotation + translation_vector;
			pose.set_xyz( core::id::AtomID( j, seq_pos ), working );
		}

	}

}

void
rigid_body_move(
	numeric::xyzMatrix<core::Real> rotation,
	numeric::xyzVector<core::Real> const & translation_vector,
	core::Real translation_scalar,
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	numeric::xyzVector<core::Real> center_of_rotation
	/* = numeric::xyzVector<core::Real>(std::numeric_limits<double>::quiet_NaN(), 0, 0) */) {

	numeric::xyzVector<core::Real> translate = translation_vector;
	translate *= translation_scalar;

	rigid_body_move( rotation, translate, pose, subset, center_of_rotation );
}



} // namespace pose_manipulation
} //toolbox
} //protocols
