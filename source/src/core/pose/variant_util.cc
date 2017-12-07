// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/variant_util.cc
/// @brief  Pose class utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Rhiju Das, Steven Lewis, Vikram K. Mulligan


// Unit header
#include <core/pose/variant_util.hh>
#include <core/pose/extra_pose_info_util.hh>

// Package headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/PositionConservedResiduesStore.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/carbohydrates/util.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/util.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/Exceptions.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.string.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/string_constants.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <boost/functional/hash.hpp>

namespace core {
namespace pose {

static basic::Tracer TR( "core.pose.variant_util" );

core::conformation::ResidueOP
remove_variant_type_from_residue(
	core::conformation::Residue const & old_rsd,
	core::chemical::VariantType const variant_type,
	pose::Pose const & pose )
{
	if ( !old_rsd.has_variant_type( variant_type ) ) return old_rsd.clone();

	// the type of the desired variant residue
	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( old_rsd.type().mode() ) );
	core::chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_removed( old_rsd.type(), variant_type ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( new_rsd_type, old_rsd, pose.conformation() ) );
	core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	if ( old_rsd.nchi() == new_rsd_type.nchi() ) {
		for ( Size chino=1; chino <= old_rsd.nchi(); ++chino ) {
			new_rsd->set_chi( chino, old_rsd.chi( chino ) );
		}
	} else {
		TR << "The chi angles will not be updated and your dunbrack score for this rotamer will be huge; "
			"this function is only meant to add a variant type to a residue of the same type" << std::endl;
	}

	return new_rsd;
}


conformation::ResidueOP
add_variant_type_to_residue(
	conformation::Residue const & old_rsd,
	chemical::VariantType const variant_type,
	pose::Pose const & pose )
{
	if ( old_rsd.has_variant_type( variant_type ) ) return old_rsd.clone();

	// the type of the desired variant residue
	chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( old_rsd.type().mode() ) );
	chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_added( old_rsd.type(), variant_type ) );
	conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( new_rsd_type, old_rsd, pose.conformation() ) );
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	if ( old_rsd.nchi() == new_rsd_type.nchi() ) {
		for ( Size chino=1; chino <= old_rsd.nchi(); ++chino ) {
			new_rsd->set_chi( chino, old_rsd.chi( chino ) );
		}
	} else {
		TR << "The chi angles will not be updated and your dunbrack score for this rotamer will be huge; "
			"this function is only meant to add a variant type to a residue of the same type" << std::endl;
	}

	return new_rsd;
}


/// @details E.g., make a terminus variant, and replace the original in pose.
/// @note This copies any atoms in common between old and new residues, rebuilding the others.
void
add_variant_type_to_pose_residue(
	pose::Pose & pose,
	chemical::VariantType const variant_type,
	Size const seqpos )
{
	runtime_assert( seqpos != 0 );
	if ( pose.residue( seqpos ).has_variant_type( variant_type ) ) return;

	conformation::Residue const & old_rsd( pose.residue( seqpos ) );

	// the type of the desired variant residue
	chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( pose.residue_type( seqpos ).mode() ) );
	chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_added( old_rsd.type(), variant_type ) );

	core::pose::replace_pose_residue_copying_existing_coordinates( pose, seqpos, new_rsd_type );

	// update connections
	for ( Size i_con=1; i_con<=pose.conformation().residue_type(seqpos).n_possible_residue_connections(); ++i_con ) {
		if ( pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con) != 0 ) {
			Size connected_seqpos = pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con);
			Size connected_id = pose.residue(seqpos).connect_map(i_con).connid();
			pose.conformation().update_noncanonical_connection(seqpos, i_con, connected_seqpos, connected_id);
		}
	}

	if ( variant_type == core::chemical::CUTPOINT_LOWER || variant_type == core::chemical::CUTPOINT_UPPER ) {
		update_cutpoint_virtual_atoms_if_connected( pose, seqpos, true );
	}
}


/// @details E.g., remove a terminus variant, and replace the original in pose.
/// @note This copies any atoms in common between old and new residues, rebuilding the others.
void
remove_variant_type_from_pose_residue(
	pose::Pose & pose,
	chemical::VariantType const variant_type,
	Size const seqpos )
{
	if ( !pose.residue( seqpos ).has_variant_type( variant_type ) ) return;

	conformation::Residue const & old_rsd( pose.residue( seqpos ) );

	// the type of the desired variant residue
	chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( pose.residue_type( seqpos ).mode() ) );
	chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_removed( old_rsd.type(), variant_type ) );

	core::pose::replace_pose_residue_copying_existing_coordinates( pose, seqpos, new_rsd_type );

	// update connections
	for ( Size i_con=1; i_con<=pose.conformation().residue_type(seqpos).n_possible_residue_connections(); ++i_con ) {
		if ( pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con) != 0 ) {
			Size connected_seqpos = pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con);
			Size connected_id = pose.residue(seqpos).connect_map(i_con).connid();
			pose.conformation().update_noncanonical_connection(seqpos, i_con, connected_seqpos, connected_id);
		}
	}
}

void
add_lower_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
) {
	add_variant_type_to_pose_residue( pose, chemical::LOWER_TERMINUS_VARIANT, seqpos );
}

void
add_upper_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
) {
	add_variant_type_to_pose_residue( pose, chemical::UPPER_TERMINUS_VARIANT, seqpos );
}

void
remove_lower_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
) {
	core::pose::remove_variant_type_from_pose_residue( pose, chemical::LOWER_TERMINUS_VARIANT, seqpos );
}

void
remove_upper_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
) {
	core::pose::remove_variant_type_from_pose_residue( pose, chemical::UPPER_TERMINUS_VARIANT, seqpos );
}

/// @brief Add cutpoint variants to all residues annotated as cutpoints in the FoldTree in the Pose.
void
correctly_add_cutpoint_variants( core::pose::Pose & pose ) {
	core::kinematics::FoldTree const & tree( pose.fold_tree() );
	for ( core::Size i = 1; i < pose.size(); ++i ) { // Less than because cutpoints are between i and i+1
		if ( tree.is_cutpoint( i ) ) {
			correctly_add_cutpoint_variants( pose, i, false );
		}
	}
}

// AMW TODO assume no foldtree perturbation
void
correctly_add_2prime_connection_variants( pose::Pose & pose, Size const twoprime_res, Size const next_res ) {

	// twoprime_res can't have any other 2prime variants
	// [none at the mo']
	using namespace core::chemical;
	using namespace core::id;

	correctly_remove_variants_incompatible_with_upper_cutpoint_variant( pose, next_res );

	// AMW: positioning TWO PRIME 'cutpoint phosphate torsions'
	// since this is all about OP1/OP2, this actually works fine... I think. It might run into issues because of how
	// it prepends a residue. yeah, major issues w/ the subsequent branch_conn scoring.

	//if ( pose.residue_type( twoprime_res ).is_RNA() )  rna::position_cutpoint_phosphate_torsions( pose, twoprime_res, next_res );

	// Manually reposition OP2, OP1 on next_res.
	AtomID aidOP1( pose.residue( next_res ).atom_index("OP1"),  next_res );
	AtomID aidOP2( pose.residue( next_res ).atom_index("OP2"),  next_res );
	AtomID aidP( pose.residue( next_res ).atom_index("P"),  next_res );
	AtomID aidO5P( pose.residue( next_res ).atom_index("O5'"),  next_res );
	AtomID aidC5P( pose.residue( next_res ).atom_index("C5'"),  next_res );

	Vector const & O2P_xyz( pose.residue( twoprime_res ).xyz( "O2'" ) );
	Vector LOWER_xyz( pose.residue( next_res ).lower_connect().icoor().build( pose.residue( next_res ), pose.conformation() ) );
	Vector const & P_xyz( pose.residue( next_res ).xyz( "P" ) );
	Vector const & O5P_xyz( pose.residue( next_res ).xyz( "O5'" ) );
	Vector const & C5P_xyz( pose.residue( next_res ).xyz( "C5'" ) );

	using namespace numeric::conversions;

	Real O2P_torsion_correction = numeric::dihedral_degrees( O2P_xyz, P_xyz, O5P_xyz, C5P_xyz ) - numeric::dihedral_degrees( LOWER_xyz, P_xyz, O5P_xyz, C5P_xyz );
	Real torsion_OP2 = degrees( pose.conformation().torsion_angle( aidOP2, aidP, aidO5P, aidC5P ) ) + O2P_torsion_correction;
	pose.conformation().set_torsion_angle( aidOP2, aidP, aidO5P, aidC5P, radians( torsion_OP2 ) );

	remove_variant_type_from_pose_residue( pose, VIRTUAL_RIBOSE, twoprime_res );
	if ( !pose.residue( twoprime_res ).has_variant_type( C2_BRANCH_POINT ) ) {
		add_variant_type_to_pose_residue( pose, C2_BRANCH_POINT, twoprime_res   );
	}
	add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, next_res );

	// important -- to prevent artificial penalty from steric clash.
	// AMW: this is a horrifying inline of declare_cutpoint_chemical_bond, but for 2' to 5'.
	using namespace core::conformation;
	// Need to clear out any chemical bonds that might have been previously tied to upper/lower of these residues.
	// check simple loop, like  in get_upper_cutpoint_partner_for_lower() in chainbreak_util.hh
	Residue const & lower_rsd( pose.conformation().residue( twoprime_res ) );
	for ( Size k = 1; k <= lower_rsd.connect_map_size(); k++ ) {
		if ( lower_rsd.residue_connect_atom_index( k ) != lower_rsd.atom_index("O2'") ) continue;
		Size upper( lower_rsd.connected_residue_at_resconn( k ) );
		if ( upper == 0 ) continue;
		Residue const & upper_rsd( pose.conformation().residue( upper ) ); // upper residue.
		Size const m = lower_rsd.residue_connection_conn_id( k );
		runtime_assert( upper_rsd.residue_connect_atom_index( m ) == upper_rsd.lower_connect_atom() );
		runtime_assert( upper_rsd.connected_residue_at_resconn( m ) == twoprime_res );
		//upper_rsd.mark_connect_incomplete( m );
		//lower_rsd.mark_connect_incomplete( k );
		pose.conformation().sever_chemical_bond( twoprime_res, k, upper, m );
	}

	Residue const & upper_rsd( pose.conformation().residue( next_res ) );
	for ( Size k = 1; k <= upper_rsd.connect_map_size(); k++ ) {
		if ( upper_rsd.residue_connect_atom_index( k ) != upper_rsd.lower_connect_atom() ) continue;
		Size lower( upper_rsd.connected_residue_at_resconn( k ) );
		if ( lower == 0 ) continue;
		Residue const & lower_rsd( pose.conformation().residue( lower ) ); // lower residue.
		Size const m = upper_rsd.residue_connection_conn_id( k );
		runtime_assert( lower_rsd.residue_connect_atom_index( m ) == lower_rsd.upper_connect_atom() );
		runtime_assert( lower_rsd.connected_residue_at_resconn( m ) == next_res );
		//lower_rsd.mark_connect_incomplete( m );
		//upper_rsd.mark_connect_incomplete( k );
		pose.conformation().sever_chemical_bond( next_res, k, lower, m );
	}

	TR << "AMW NOTE: O2' res " << lower_rsd.name() << " connect map is " << lower_rsd.connect_map_size() << std::endl;
	pose.conformation().declare_chemical_bond(
		twoprime_res,
		"O2'",
		next_res,
		upper_rsd.atom_name( upper_rsd.lower_connect_atom() ) );

	TR << "AMW NOTE: O2' res " << lower_rsd.name() << " connect map is " << lower_rsd.connect_map_size() << std::endl;
}

/// @brief Add CUTPOINT_LOWER and CUTPOINT_UPPER types to two residues, remove incompatible types, and declare
/// a chemical bond between them.
/// @param[in,out] pose The pose to modify.
/// @param[in] cutpoint_res The index of the CUTPOINT_LOWER residue.
/// @param[in] check_fold_tree If true, a check is performed to confirm that the residues in question represent a
/// cutpoint in the foldtree in the pose.
/// @param[in] next_res_in The index of the CUTPOINT_UPPER residue.  If not provided, or if set to 0, this defaults
/// to the cutpoint_res + 1 residue.  Must be specified for cyclic geometry.
void
correctly_add_cutpoint_variants(
	core::pose::Pose & pose,
	Size const cutpoint_res,
	bool const check_fold_tree /* = true*/,
	Size const next_res_in /* = 0 */
) {
	using namespace core::chemical;

	Size const next_res = ( next_res_in == 0 ) ? ( cutpoint_res + 1 ) : next_res_in; // user might specify a different "next_res" to cyclize.

	if ( check_fold_tree ) runtime_assert( pose.fold_tree().is_cutpoint( cutpoint_res ) );

	correctly_remove_variants_incompatible_with_lower_cutpoint_variant( pose, cutpoint_res );
	correctly_remove_variants_incompatible_with_upper_cutpoint_variant( pose, next_res );

	if ( pose.residue_type( cutpoint_res ).is_RNA() ) rna::position_cutpoint_phosphate_torsions( pose, cutpoint_res, next_res );

	add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint_res );
	add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, next_res );

	// important -- to prevent artificial penalty from steric clash.
	declare_cutpoint_chemical_bond( pose, cutpoint_res, next_res );

	// Update positions of virtual atoms:
	if ( !pose.residue( cutpoint_res ).is_RNA() ) update_cutpoint_virtual_atoms_if_connected( pose, cutpoint_res, false );
	if ( !pose.residue( next_res     ).is_RNA() ) update_cutpoint_virtual_atoms_if_connected( pose, next_res, false );
}

/// @brief Remove variant types incompatible with CUTPOINT_LOWER from a position in a pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
/// @param[in,out] pose The pose on which to operate.
/// @param[in] res_index The index of the residue on which to operate.
void
correctly_remove_variants_incompatible_with_lower_cutpoint_variant(
	core::pose::Pose & pose,
	Size const res_index
) {
	remove_variant_type_from_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, res_index );
	remove_variant_type_from_pose_residue( pose, core::chemical::THREE_PRIME_PHOSPHATE, res_index );
	remove_variant_type_from_pose_residue( pose, core::chemical::C_METHYLAMIDATION, res_index );
}

/// @brief Remove variant types incompatible with CUTPOINT_UPPER from a position in a pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
/// @param[in,out] pose The pose on which to operate.
/// @param[in] res_index The index of the residue on which to operate.
void
correctly_remove_variants_incompatible_with_upper_cutpoint_variant(
	core::pose::Pose & pose,
	Size const res_index
) {
	remove_variant_type_from_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, res_index );
	remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, res_index );
	remove_variant_type_from_pose_residue( pose, core::chemical::FIVE_PRIME_PHOSPHATE, res_index );
	remove_variant_type_from_pose_residue( pose, core::chemical::N_ACETYLATION, res_index);
}

/// @brief returns true if the given residue in the pose is a chain ending or has upper/lower terminal variants
bool
pose_residue_is_terminal( Pose const & pose, Size const resid )
{
	return ( is_lower_terminus( pose, resid ) || is_upper_terminus( pose, resid ) );
}

/// @brief checks to see if this is a lower chain ending more intelligently than just checking residue variants
bool
is_lower_terminus( pose::Pose const & pose, Size const resid )
{
	return ( ( resid == 1 ) || //loop starts at first residue
		( ! pose.residue(resid).is_polymer() ) || // this residue isn't a polymer
		( ! pose.residue( resid-1 ).is_protein() ) || //residue before start is not protein
		( pose.chain( resid-1 ) != pose.chain( resid ) ) || // residues before start are on another chain
		( pose.residue( resid ).is_lower_terminus() ) ); // start of residue is lower terminus
}

/// @brief checks to see if this is a lower chain ending more intelligently than just checking residue variants
bool
is_upper_terminus( pose::Pose const & pose, Size const resid )
{
	return ( ( resid == pose.size() ) || // loop end at last residue
		( !pose.residue( resid ).is_polymer() ) || // this residue isn't a polymer
		( !pose.residue( resid+1 ).is_protein() ) || // residue after end is not protein
		( pose.chain( resid+1 ) != pose.chain( resid ) ) || // residues before start is other chain
		( pose.residue( resid ).is_upper_terminus() ) ); // explicit terminus variant @ end of loop
}

////////////////////////////////////////////////////////////////////
void
fix_up_residue_type_variants_at_strand_end( pose::Pose & pose, Size const res ) {

	using namespace core::chemical;
	using namespace core::pose::full_model_info;

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > const chains_full = get_chains_full( pose );

	// Could this be a chainbreak (cutpoint_closed )?
	TR.Debug << "checking for cutpoint after append: " << res << " " << res_list[ res ]  << " " << cutpoint_open_in_full_model.size() << std::endl;

	if ( res < pose.size() &&
			res_list[ res ] + 1 == res_list[ res + 1 ] &&
			! cutpoint_open_in_full_model.has_value( res_list[ res ]) ) {

		if ( pose.residue_type( res ).has_variant_type( CUTPOINT_LOWER ) &&
				pose.residue_type( res + 1 ).has_variant_type( CUTPOINT_UPPER ) ) return;

		// can happen after additions
		core::pose::correctly_add_cutpoint_variants( pose, res );

		// leave virtual riboses in this should actually get instantiated by the modeler
		// remove_variant_type_from_pose_residue( pose, chemical::VIRTUAL_RIBOSE, res );

	} else {

		// can happen after deletions
		remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, res );

		// proteins...
		if ( pose.residue_type( res ).is_protein() ) {
			if ( res_list[ res ] < full_model_info.size() &&
					chains_full[ res_list[ res ] + 1 ] == chains_full[ res_list[ res ] ] &&
					( res == pose.size()  || res_list[ res ] + 1 < res_list[ res + 1 ] ) ) {
				remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS_VARIANT, res );
				add_variant_type_to_pose_residue( pose, C_METHYLAMIDATION, res );
			} else {
				remove_variant_type_from_pose_residue( pose, C_METHYLAMIDATION, res );
				add_variant_type_to_pose_residue( pose, UPPER_TERMINUS_VARIANT, res );
			}
		}

	}

}

////////////////////////////////////////////////////////////////////
void
fix_up_residue_type_variants_at_strand_beginning( pose::Pose & pose, Size const res ) {

	using namespace core::chemical;
	using namespace core::pose::full_model_info;

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > const chains_full = get_chains_full( pose );

	// Could this be a chainbreak (cutpoint_closed )?

	TR.Debug << "checking for cutpoint after prepend: " << res << " " << res_list[ res ] << " " << cutpoint_open_in_full_model.size() << std::endl;

	if ( res > 1 &&
			res_list[ res ] - 1 == res_list[ res - 1 ] &&
			! cutpoint_open_in_full_model.has_value( res_list[ res - 1 ])  ) {

		if ( pose.residue_type( res - 1 ).has_variant_type( CUTPOINT_LOWER ) &&
				pose.residue_type( res     ).has_variant_type( CUTPOINT_UPPER ) ) return;

		// can happen after additions
		core::pose::correctly_add_cutpoint_variants( pose, res - 1 );
	} else {
		// can happen after additions
		if ( pose.residue_type( res ).is_RNA() &&
				!pose.residue_type( res ).has_variant_type( FIVE_PRIME_PHOSPHATE ) ) {
			add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, res );
		}

		// can happen after deletions
		remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, res );

		// proteins...
		if ( pose.residue_type( res ).is_protein() ) {
			if ( res_list[ res ] > 1 &&
					chains_full[ res_list[ res ] - 1 ] == chains_full[ res_list[ res ] ] &&
					( res == 1 || res_list[ res ] - 1 > res_list[ res - 1 ] ) ) {
				remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS_VARIANT, res );
				add_variant_type_to_pose_residue( pose, N_ACETYLATION, res );
			} else {
				remove_variant_type_from_pose_residue( pose, N_ACETYLATION, res );
				add_variant_type_to_pose_residue( pose, LOWER_TERMINUS_VARIANT, res );
			}
		}
	}
}

////////////////////////////////////////////////////////////////////
void
fix_up_residue_type_variants_at_floating_base( pose::Pose & pose, Size const res ) {

	using namespace full_model_info;

	if ( !pose.residue_type(res ).is_RNA() ) return;
	remove_variant_type_from_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, res );
	remove_variant_type_from_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, res );


	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > const & sample_res = full_model_info.sample_res();
	utility::vector1< Size > const & sample_sugar_res = full_model_info.rna_sample_sugar_res();

	if ( !sample_res.has_value( res_list[ res ] ) &&
			!sample_sugar_res.has_value( res_list[ res ] ) ) return;

	if ( res > 1 &&
			res_list[ res ] - 1 == res_list[ res - 1 ] &&
			! cutpoint_open_in_full_model.has_value( res_list[ res - 1 ])  ) return;

	if ( res < pose.size() &&
			res_list[ res ] + 1 == res_list[ res + 1 ] &&
			! cutpoint_open_in_full_model.has_value( res_list[ res ]) ) return;

	if ( pose.residue_type( res ).has_variant_type( core::chemical::FIVE_PRIME_PHOSPHATE ) )  return;
	if ( pose.residue_type( res ).has_variant_type( core::chemical::THREE_PRIME_PHOSPHATE ) ) return;

	add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, res );

}

////////////////////////////////////////////////////////////////////
void
update_block_stack_variants( pose::Pose & pose, Size const & n ) {
	using namespace core::chemical;
	using namespace core::pose::full_model_info;
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & block_stack_above_res = full_model_info.rna_block_stack_above_res();
	utility::vector1< Size > const & block_stack_below_res = full_model_info.rna_block_stack_below_res();

	if ( block_stack_above_res.has_value( res_list[ n ] ) ) {
		add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, n );
	} else {
		runtime_assert( !pose.residue_type( n ).has_variant_type( BLOCK_STACK_ABOVE ) );
	}
	if ( block_stack_below_res.has_value( res_list[ n ] ) ) {
		add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, n );
	} else {
		runtime_assert( !pose.residue_type( n ).has_variant_type( BLOCK_STACK_BELOW ) );
	}

}

////////////////////////////////////////////////////////////////////
void
fix_up_residue_type_variants(
#ifndef GL_GRAPHICS
	pose::Pose & pose
#else
	pose::Pose & pose_to_fix 
#endif
) {

	using namespace core::chemical;
	using namespace core::pose::full_model_info;
#ifdef GL_GRAPHICS
	pose::Pose pose = pose_to_fix; // costly, but prevents seg fault with graphics.
#endif

	auto const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = full_model_info.res_list();

	for ( Size n = 1; n <= pose.size(); n++ ) {

		// Are we at a strand beginning?
		bool const at_strand_beginning = ( n == 1 || pose.fold_tree().is_cutpoint( n-1 ) );
		if ( at_strand_beginning ) {
			fix_up_residue_type_variants_at_strand_beginning( pose, n );
		} else { // make sure there is nothing crazy here
			remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS_VARIANT, n );
			remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, n );
			remove_variant_type_from_pose_residue( pose, FIVE_PRIME_PHOSPHATE, n );
			remove_variant_type_from_pose_residue( pose, VIRTUAL_RIBOSE, n );
			runtime_assert( !pose.residue_type( n ).has_variant_type( CUTPOINT_UPPER ) );
		}

		// Look for strand ends.
		bool const at_strand_end = ( n == pose.size() || pose.fold_tree().is_cutpoint( n ) );
		if ( at_strand_end ) {
			fix_up_residue_type_variants_at_strand_end( pose, n );
		} else {
			remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS_VARIANT, n );
			remove_variant_type_from_pose_residue( pose, VIRTUAL_RIBOSE, n );
			remove_variant_type_from_pose_residue( pose, THREE_PRIME_PHOSPHATE, n );
			runtime_assert( !pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) );
		}

		// check for floating_base
		if ( at_strand_end && at_strand_beginning ) fix_up_residue_type_variants_at_floating_base( pose,  n );

		update_block_stack_variants( pose, n );
	}

	for ( auto const & cyclize_pair : full_model_info.cyclize_res() ) {
		Size const first = cyclize_pair.first;
		Size const second = cyclize_pair.second;
		// Bigger number first... actually use max/min in case wrong order.
		if ( res_list.contains( first ) && res_list.contains( second ) ) {
			TR << "OK, now gonna add cutpoint variants (correctly!)" << std::endl;
			core::pose::correctly_add_cutpoint_variants( pose, std::max( res_list.index( first ), res_list.index( second ) ), false, std::min( res_list.index( first ), res_list.index( second ) ) );
		}
	}

	for ( auto const & twoprime_pair : full_model_info.twoprime_res() ) {
		Size const first = twoprime_pair.first;
		Size const second = twoprime_pair.second;
		// first one ==  the one with the 2prime variant
		// Ugh, but this function sorts them automatically.
		// Since TYPICALLY smaller-to-larger connections are polymeric
		// assume larger-to-smaller
		if ( res_list.contains( second ) && res_list.contains( first ) ) {
			TR << "OK, now gonna add 2prime variants (correctly!) " << res_list.index( first ) << " " << res_list.index( second ) << std::endl;
			core::pose::correctly_add_2prime_connection_variants( pose, res_list.index( second ), res_list.index( first ) );
		}
	}

#ifdef GL_GRAPHICS
	// Just copying the conformation() makes sure that other objects (such as other_pose_list) don't get cloned --
	//  can be important if external functions are holding OPs to those objects.
	pose_to_fix.conformation() = pose.conformation();
	pose_to_fix.pdb_info( pose.pdb_info() ); // silly -- ensures that PDBInfo is not flagged as 'obsolete'.
#endif
}

} // pose
} // core
