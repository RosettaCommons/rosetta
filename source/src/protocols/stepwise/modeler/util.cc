// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file stepwise/modeler/util
/// @brief Utility functions specific for stepwise modeler
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/rna/denovo/secstruct_legacy/RNA_SecStructLegacyInfo.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/copydofs/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/loop_graph/LoopGraph.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/methods/chainbreak_util.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/stream_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

using namespace core;
using namespace core::pose;
using namespace core::pose::full_model_info;

using ObjexxFCL::string_of;
using utility::operator<<;
using core::kinematics::FoldTree;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.util" );

//////////////////////////////////////////////////////////////////////////
//
// Recently cleaned up to only include critical pose cut/paste/merge/split
//  functionality for stepwise assembly/monte-carlo.
//
// Some very very useful functions in here.
//
//////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace stepwise {
namespace modeler {

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
Size
make_cut_at_moving_suite( core::kinematics::FoldTree & fold_tree, Size const & moving_suite ){

	fold_tree.new_jump( moving_suite, moving_suite + 1, moving_suite );
	return find_jump_number_at_suite( fold_tree, moving_suite );

}

///////////////////////////////////////////////////////////////////////
Size
make_cut_at_moving_suite( pose::Pose & pose, Size const & moving_suite ){

	core::kinematics::FoldTree fold_tree = pose.fold_tree();

	Size jump_number( 0 );
	if ( fold_tree.is_cutpoint( moving_suite ) ) { // already a cutpoint there
		return find_jump_number_at_suite( fold_tree, moving_suite );
	} else {
		jump_number = make_cut_at_moving_suite( fold_tree, moving_suite );
	}

	pose.fold_tree( fold_tree );

	return jump_number;
}

///////////////////////////////////////////////////////////////////////
Size
find_jump_number_at_suite( kinematics::FoldTree const & fold_tree, Size const & moving_suite ){

	int const i( moving_suite ), j( moving_suite + 1 );
	for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
		if ( fold_tree.upstream_jump_residue( n ) == i && fold_tree.downstream_jump_residue( n ) == j ) return n;
		if ( fold_tree.upstream_jump_residue( n ) == j && fold_tree.downstream_jump_residue( n ) == i ) return n;
	}

	utility_exit_with_message( "Problem with jump number" );
	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////
Size
look_for_unique_jump_to_moving_res( kinematics::FoldTree const & fold_tree, Size const & i ){

	for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
		if ( fold_tree.downstream_jump_residue( n ) == int( i )  ) {
			return n;
		}
	}


	Size num_jump( 0 ), jump_idx( 0 );
	for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
		if ( fold_tree.upstream_jump_residue( n ) == int( i ) || fold_tree.downstream_jump_residue( n ) == int( i ) ) {
			jump_idx = n;
			num_jump++;
		}
	}
	runtime_assert( num_jump == 1 );
	return jump_idx;
}

///////////////////////////////////////////////////////////////////////////////////////
bool
is_cutpoint_closed( pose::Pose const & pose, Size const seq_num ){
	runtime_assert( seq_num > 0 );
	runtime_assert( seq_num <= pose.size() );
	if ( pose.residue_type( seq_num  ).has_variant_type( chemical::CUTPOINT_LOWER )  ) {
		Size const j = core::scoring::methods::get_upper_cutpoint_partner_for_lower( pose, seq_num );
		if ( !pose.residue_type( j ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
			utility_exit_with_message( "found CUTPOINT_LOWER with no CUTPOINT_UPPER" );
		}
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_cutpoint_closed( pose::Pose const & pose ){
	utility::vector1< Size > cutpoint_closed;
	for ( Size seq_num = 1; seq_num < pose.size(); seq_num++ ) {
		if ( is_cutpoint_closed( pose, seq_num ) ) cutpoint_closed.push_back( seq_num );
	}
	return cutpoint_closed;
}


///////////////////////////////////////////////////////////////////////
utility::vector1< Size >
merge_vectors( utility::vector1< Size > const & vec1,
	utility::vector1< Size > const & vec2 ){

	std::map< Size, bool > silly_map;
	for ( Size n = 1; n <= vec1.size(); n++ ) silly_map[ vec1[n] ] = true;
	for ( Size n = 1; n <= vec2.size(); n++ ) silly_map[ vec2[n] ] = true;

	utility::vector1< Size > merged_vec;
	for ( std::map<Size,bool>::iterator it = silly_map.begin(), end = silly_map.end(); it != end; ++it ) {
		merged_vec.push_back( it->first );
	}
	return merged_vec;

}

////////////////////////////////////////////////////////////////////////////////////////////////
void
get_jump_partners_from_pose( utility::vector1< Size > & jump_partners1,
	utility::vector1< Size > & jump_partners2,
	utility::vector1< std::string > & jump_atoms1,
	utility::vector1< std::string > & jump_atoms2,
	pose::Pose const & pose,
	utility::vector1< Size > const & working_res ){

	using namespace core::kinematics;

	FoldTree const & f = pose.fold_tree();
	for ( Size n = 1; n <= f.num_jump(); n++ ) {
		Size const j1 = f.upstream_jump_residue( n );
		Size const j2 = f.downstream_jump_residue( n );
		jump_partners1.push_back(  working_res[ j1 ] );
		jump_partners2.push_back(  working_res[ j2 ] );
		jump_atoms1.push_back( f.upstream_atom( n ) );
		jump_atoms2.push_back( f.downstream_atom( n ) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
void
get_endpoints_from_pose( utility::vector1< Size > & endpoints,
	pose::Pose const & pose,
	utility::vector1< Size > const & working_res ){

	for ( Size i = 1; i < pose.size(); i++ ) {
		if ( pose.fold_tree().is_cutpoint( i ) ) {
			endpoints.push_back( working_res[ i ] );
		}
	}
	endpoints.push_back( working_res[ pose.size() ] );

}


////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
map_to_local_numbering( utility::vector1< Size > const & vec,
	utility::vector1< Size > const & working_res ){
	utility::vector1< Size > vec_new;
	for ( Size n = 1; n <= vec.size(); n++ ) {
		vec_new.push_back( working_res.index( vec[n] ) );
	}
	return vec_new;
}

////////////////////////////////////////////////////////////////////////////////////////////////
std::map< Size, Size >
get_res_map( utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & source_res ){

	std::map< Size, Size > res_map;

	for ( Size n = 1; n <= working_res.size(); n++ ) {
		if ( !source_res.has_value( working_res[n] ) ) continue;
		res_map[ n ] = source_res.index( working_res[n] );
	}

	return res_map;
}


////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
merge_disjoint_vectors( utility::vector1< Size > const & res_vector1,
	utility::vector1< Size > const & res_vector2 ){

	utility::vector1< Size > res_vector = res_vector1;

	for ( Size n = 1; n <= res_vector2.size(); n++ ) {
		runtime_assert( ! res_vector.has_value( res_vector2[n] ) );
		res_vector.push_back( res_vector2[ n ] );
	}

	std::sort( res_vector.begin(), res_vector.end() ); // hope this works.

	return res_vector;
}

////////////////////////////////////////////////////////////////////////////////////////////////
bool
find_root_without_virtual_ribose( kinematics::FoldTree const & f, pose::Pose const & pose ){
	for ( Size n = 1; n <= f.nres(); n++ ) {
		if ( f.possible_root( n ) && !pose.residue_type( n ).has_variant_type( chemical::VIRTUAL_RIBOSE ) ) {
			return n;
		}
	}
	utility_exit_with_message( "Fail: find_root_without_virtual_ribose" );
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////
Size
find_first_root_residue( kinematics::FoldTree const & f,
	utility::vector1< Size > const & working_res_subset,
	utility::vector1< Size > const & working_res ){
	runtime_assert( f.nres() == working_res.size() );
	bool found_root( false );
	for ( Size n = 1; n <= f.nres(); n++ ) {
		if ( working_res_subset.has_value( working_res[ n ] ) &&
				f.possible_root( n ) ) {
			found_root = true;
			return n;
		}
	}
	runtime_assert( found_root );
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void
reroot_based_on_full_model_info( pose::Pose & pose ) {
	utility::vector1< Size > root_partition_res;
	for ( Size n = 1; n <= pose.size(); n++ ) root_partition_res.push_back( n );
	reroot_based_on_full_model_info( pose, root_partition_res );
}

////////////////////////////////////////////////////////////////////////////////////////////////
void
reroot_based_on_full_model_info( pose::Pose & pose,
	utility::vector1< Size > const & root_partition_res ){
	// first look for any user-defined roots.
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & preferred_root_res = full_model_info.preferred_root_res();
	utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > const & working_res = full_model_info.working_res();
	core::pose::reroot( pose, root_partition_res, res_list, preferred_root_res, fixed_domain_map,
		cutpoint_open_in_full_model, working_res );
}



////////////////////////////////////////////////////////////////////////////////////////////////
void
merge_in_other_pose_by_jump( pose::Pose & pose, pose::Pose const & pose2,
	Size const lower_merge_res, Size const upper_merge_res  ){
	runtime_assert( lower_merge_res < upper_merge_res );
	merge_in_other_pose( pose, pose2,
		lower_merge_res, upper_merge_res,
		false /*connect_residues_by_bond*/ );
}

////////////////////////////////////////////////////////////////////////////////////////////////
void
merge_in_other_pose_by_bond( pose::Pose & pose, pose::Pose const & pose2, Size const merge_res ){
	Size const lower_merge_res( merge_res ), upper_merge_res( merge_res + 1 );
	merge_in_other_pose( pose, pose2,
		lower_merge_res, upper_merge_res,
		true /*connect_residues_by_bond*/ );
}

////////////////////////////////////////////////////////////////////////////////////////////////
void
merge_in_other_pose( pose::Pose & pose, pose::Pose const & pose2,
	Size const lower_merge_res, Size const upper_merge_res,
	bool const connect_residues_by_bond ) {

	using namespace core::pose::datacache;
	using namespace core::pose::full_model_info;

	if ( connect_residues_by_bond ) runtime_assert( upper_merge_res == lower_merge_res + 1 );

	FullModelInfo & full_model_info = nonconst_full_model_info( pose );

	core::pose::Pose pose_scratch;
	utility::vector1< Size > const new_res_list = merge_two_poses_using_full_model_info( pose_scratch, pose, pose2,
		lower_merge_res, upper_merge_res,
		connect_residues_by_bond );
	pose.conformation() = pose_scratch.conformation();

	full_model_info.set_res_list( new_res_list );
	full_model_info.add_submotif_info( const_full_model_info( pose2 ).submotif_info_list() );
	full_model_info.update_submotif_info_list();
	update_pose_objects_from_full_model_info( pose ); // for output pdb or silent file -- residue numbering.

}

////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
merge_two_poses_using_full_model_info( pose::Pose & pose,
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	Size const lower_merge_res,
	Size const upper_merge_res,
	bool const connect_residues_by_bond ) {

	using namespace core::pose::full_model_info;
	using namespace core::pose::datacache;
	using namespace basic::datacache;

	// get working_residue information from each pose
	utility::vector1< Size > const & working_res1 = get_res_list_from_full_model_info_const( pose1 );
	utility::vector1< Size > const & working_res2 = get_res_list_from_full_model_info_const( pose2 );

	if ( connect_residues_by_bond ) runtime_assert( upper_merge_res == lower_merge_res + 1 );

	return merge_two_poses( pose, pose1, pose2, working_res1, working_res2,
		lower_merge_res, upper_merge_res, connect_residues_by_bond );

}


////////////////////////////////////////////////////////////////////////////////////////////////
// Note that following is really general, and potentially very valuable, for all stepwise approaches
//  and even in pose setup for fragment assembly of RNA or proteins.
// Consider including in a GeneralPoseSetup class.
//
// Trying to get correct fold-tree handling, and variants.
// Also trying to properly handle sequence reorderings (which can get pretty complicated in general )
//
utility::vector1< Size >
merge_two_poses( pose::Pose & pose,
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	utility::vector1< Size > const & working_res1,
	utility::vector1< Size > const & working_res2,
	Size const lower_merge_res,
	Size const upper_merge_res,
	bool const connect_residues_by_bond,
	bool const fix_first_pose /* = true */) {

	using namespace core::kinematics;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pose::copydofs;

	if ( connect_residues_by_bond ) runtime_assert( upper_merge_res == lower_merge_res + 1 );

	if ( working_res1.has_value( lower_merge_res ) ) {
		runtime_assert( working_res2.has_value( upper_merge_res ) );
	} else {
		runtime_assert( working_res2.has_value( lower_merge_res   ) );
		runtime_assert( working_res1.has_value( upper_merge_res ) );
		TR.Debug << "merge_two_poses: order is switched.  " << std::endl;
		// order is switched. No problem. Well... this could probably be avoided if we instead
		// specify merge_res1 and merge_res2, and allow them to be in arbitrary sequential order.
		return merge_two_poses( pose, pose2, pose1, working_res2, working_res1,
			lower_merge_res, upper_merge_res, connect_residues_by_bond,
			! fix_first_pose );
	}

	// define full working res -- union of the working res of the individual poses.
	utility::vector1< Size > const working_res = merge_disjoint_vectors( working_res1, working_res2 );

	for ( Size n = 1; n <= working_res.size(); n++ ) {
		Size const k = working_res[ n ];
		ResidueOP rsd;
		bool after_cutpoint( false );
		bool rsd_from_2(false);
		if ( working_res1.has_value( k ) ) {
			Size const j = working_res1.index( k );
			if ( j == 1 || pose1.fold_tree().is_cutpoint( j-1 ) ) after_cutpoint = true;
			rsd = pose1.residue( j ).clone();
			rsd_from_2=false;
		} else {
			runtime_assert( working_res2.has_value( k ) );
			Size const j = working_res2.index( k );
			if ( j == 1 || pose2.fold_tree().is_cutpoint( j-1 ) ) after_cutpoint = true;
			rsd = pose2.residue( j ).clone();
			rsd_from_2=true;
		}

		if ( connect_residues_by_bond && k == lower_merge_res ) {
			rsd = remove_variant_type_from_residue( *rsd, chemical::UPPER_TERMINUS_VARIANT, (rsd_from_2 ? pose2 : pose1) ); // got to be safe.
			rsd = remove_variant_type_from_residue( *rsd, chemical::THREE_PRIME_PHOSPHATE, (rsd_from_2 ? pose2 : pose1) ); // got to be safe.
			rsd = remove_variant_type_from_residue( *rsd, chemical::C_METHYLAMIDATION, (rsd_from_2 ? pose2 : pose1) ); // got to be safe.
		}
		if ( connect_residues_by_bond && k == upper_merge_res ) {
			runtime_assert( after_cutpoint );
			rsd = remove_variant_type_from_residue( *rsd, chemical::LOWER_TERMINUS_VARIANT, (rsd_from_2 ? pose2 : pose1) ); // got to be safe.
			rsd = remove_variant_type_from_residue( *rsd, chemical::FIVE_PRIME_PHOSPHATE, (rsd_from_2 ? pose2 : pose1) ); // got to be safe.
			rsd = remove_variant_type_from_residue( *rsd, chemical::N_ACETYLATION, (rsd_from_2 ? pose2 : pose1) ); // got to be safe.
			after_cutpoint = false; // we're merging after all.
		}
		if ( rsd->is_polymer() && ( n == 1 || !after_cutpoint ) ) {
			pose.append_residue_by_bond( *rsd, true /* build_ideal_geometry */ );
		} else {
			pose.append_residue_by_jump( *rsd, pose.size() );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	// figure out fold tree for this merged pose --> move this to its own function.

	// figure out jumps
	utility::vector1< Size > jump_partners1, jump_partners2, endpoints, cuts;
	utility::vector1< std::string > jump_atoms1, jump_atoms2;
	get_jump_partners_from_pose( jump_partners1, jump_partners2, jump_atoms1, jump_atoms2, pose1, working_res1 );
	get_jump_partners_from_pose( jump_partners1, jump_partners2, jump_atoms1, jump_atoms2, pose2, working_res2 );
	if ( !connect_residues_by_bond ) {
		jump_partners1.push_back( lower_merge_res );
		jump_partners2.push_back( upper_merge_res );
		jump_atoms1.push_back( chemical::rna::default_jump_atom( pose1.residue_type( working_res1.index( lower_merge_res ) ) ) ); // rna specific?
		jump_atoms2.push_back( chemical::rna::default_jump_atom( pose2.residue_type( working_res2.index( upper_merge_res ) ) ) ); // rna specific?
	}
	runtime_assert( jump_partners1.size() == jump_partners2.size() );

	// figure out cuts
	get_endpoints_from_pose( endpoints, pose1, working_res1 );
	get_endpoints_from_pose( endpoints, pose2, working_res2 );
	if ( connect_residues_by_bond ) runtime_assert( endpoints.has_value( lower_merge_res ) );
	Size const last_res = working_res[ working_res.size() ];
	runtime_assert( endpoints.has_value( last_res ) );
	for ( Size n = 1; n <= endpoints.size(); n++ ) {
		if ( ( !connect_residues_by_bond || endpoints[ n ] != lower_merge_res ) &&
				endpoints[ n ] != last_res ) cuts.push_back( endpoints[ n ] );
	}

	Size num_cuts = cuts.size();
	runtime_assert( num_cuts == jump_partners1.size() );

	jump_partners1 = map_to_local_numbering( jump_partners1, working_res );
	jump_partners2 = map_to_local_numbering( jump_partners2, working_res );
	cuts = map_to_local_numbering( cuts, working_res );

	FoldTree f = core::import_pose::get_tree( pose.size(), cuts, jump_partners1, jump_partners2, jump_atoms1, jump_atoms2 );

	Size root( 0 );
	if ( fix_first_pose ) {
		root = find_first_root_residue( f, working_res1, working_res );
	} else {
		root = find_first_root_residue( f, working_res2, working_res );
	}
	runtime_assert( f.reorder( root ) );

	pose.fold_tree( f );

	// map (internal) coordinates from separate poses into merged one.
	// this is potentially dangerous, if we accumulate floating point errors, and scores change.
	std::map< Size, Size > res_map1 = get_res_map( working_res, working_res1 );
	std::map< Size, Size > res_map2 = get_res_map( working_res, working_res2 );

	copy_dofs_match_atom_names( pose, pose1, res_map1, false /*backbone_only*/, false /*side_chain_only*/, false /*ignore_virtual*/ );
	copy_dofs_match_atom_names( pose, pose2, res_map2, false /*backbone_only*/, false /*side_chain_only*/, false /*ignore_virtual*/ );

	declare_chemical_bonds_at_cutpoints( pose, pose1, working_res, working_res1 );
	declare_chemical_bonds_at_cutpoints( pose, pose2, working_res, working_res2 );

	if ( fix_first_pose ) {
		align::superimpose_pose_legacy( pose, pose1, res_map1 );
	} else {
		align::superimpose_pose_legacy( pose, pose2, res_map2 );
	}

	return working_res;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void
declare_chemical_bonds_at_cutpoints( pose::Pose & pose,
	pose::Pose const & source_pose,
	utility::vector1< Size > const & working_res )
{
	utility::vector1< Size > source_res;
	for ( Size n = 1; n <= source_pose.size(); n++ ) source_res.push_back( n );
	declare_chemical_bonds_at_cutpoints( pose, source_pose, working_res, source_res );
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Maybe should generalize this to handle all connections, not just bonds across cutpoints --
//   would be important for, e.g., disulfides.
////////////////////////////////////////////////////////////////////////////////////////////////
void
declare_chemical_bonds_at_cutpoints( pose::Pose & pose,
	pose::Pose const & source_pose,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & source_working_res )
{
	using namespace core::chemical;
	using namespace core::scoring::methods;
	for ( Size n = 1; n <= source_pose.size(); n++ ) {
		if ( source_pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ) {
			Size next_res = get_upper_cutpoint_partner_for_lower( source_pose, n );
			if ( !working_res.has_value( source_working_res[ n ] ) ) continue;
			if ( !working_res.has_value( source_working_res[ next_res ] ) ) continue;
			declare_cutpoint_chemical_bond( pose,
				working_res.index( source_working_res[ n ] ),
				working_res.index( source_working_res[ next_res ] ) );
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
apply_numbering( utility::vector1< Size > const & res,
	utility::vector1< Size > const & numbering ) {
	utility::vector1< Size > res_renumbered;
	for ( Size n = 1; n <= res.size(); n++ ) {
		runtime_assert( res[ n ] > 0 && res[ n ] <= numbering.size() );
		res_renumbered.push_back( numbering[ res[ n ] ] );
	}
	return res_renumbered;
}

////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_other_residues( utility::vector1< Size > const & res, Size const & nres ){
	utility::vector1< Size > other_res;
	for ( Size n = 1; n <= nres; n++ ) {
		if ( !res.has_value( n ) ) other_res.push_back( n );
	}
	return other_res;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// could/should also make a version slice_out_pose_and_update_full_model_info...
void
slice_out_pose( pose::Pose & pose,
	pose::Pose & sliced_out_pose,
	utility::vector1< Size > const & residues_to_delete ) {

	using namespace core::pose;
	using namespace core::pose::full_model_info;

	// need this for later.
	FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	utility::vector1< Size > const original_res_list = full_model_info.res_list();

	// piece sliced out
	slice( sliced_out_pose, pose, residues_to_delete );

	FullModelInfoOP sliced_out_full_model_info = full_model_info.clone_info();
	sliced_out_full_model_info->set_res_list( apply_numbering( residues_to_delete, original_res_list ) );
	sliced_out_full_model_info->update_submotif_info_list();
	sliced_out_full_model_info->clear_other_pose_list();
	sliced_out_pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, sliced_out_full_model_info );
	update_pose_objects_from_full_model_info( sliced_out_pose ); // for output pdb or silent file -- residue numbering.
	protocols::scoring::set_vdw_cached_rep_screen_info_from_pose( sliced_out_pose, pose );

	// remainder piece.
	utility::vector1< Size > const residues_to_retain = get_other_residues( residues_to_delete, pose.size() );

	Pose pose_scratch;
	slice( pose_scratch, pose, residues_to_retain );
	pose.conformation() = pose_scratch.conformation();

	full_model_info.set_res_list( apply_numbering( residues_to_retain, original_res_list ) );
	full_model_info.update_submotif_info_list();
	update_pose_objects_from_full_model_info( pose ); // for output pdb or silent file -- residue numbering.
	protocols::scoring::set_vdw_cached_rep_screen_info_from_pose( pose, sliced_out_pose );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// this may falter if the sliced out pose is not contiguous. That's what the boolean return value is for.
// Could also return MovingResidueCase
void
slice( pose::Pose & sliced_out_pose,
	pose::Pose const & pose,
	utility::vector1< Size > const & slice_res ){

	using namespace core::kinematics;
	using namespace core::conformation;
	using namespace core::pose::copydofs;

	Size num_five_prime_connections( 0 ), num_three_prime_connections( 0 );
	Size num_jumps_to_previous( 0 ), num_jumps_to_next( 0 );

	sliced_out_pose.conformation().clear();

	for ( Size n = 1; n <= slice_res.size(); n++ ) {

		Size const j = slice_res[ n ];
		ResidueOP rsd = pose.residue( j ).clone();

		bool const after_cutpoint = ( j == 1 || ( pose.fold_tree().is_cutpoint( j - 1 ) ) );
		if ( after_cutpoint ) rsd = remove_variant_type_from_residue( *rsd, chemical::LOWER_TERMINUS_VARIANT, pose );
		if ( !after_cutpoint && !slice_res.has_value( j - 1 ) ) num_five_prime_connections++;

		bool const jump_to_previous = check_jump_to_previous_residue_in_chain( pose, j, slice_res );
		if ( jump_to_previous ) num_jumps_to_previous++;

		bool const before_cutpoint = ( j == pose.size() || ( pose.fold_tree().is_cutpoint( j ) ) );
		if ( before_cutpoint ) rsd = remove_variant_type_from_residue( *rsd, chemical::UPPER_TERMINUS_VARIANT, pose ); // got to be safe.
		if ( !before_cutpoint && !slice_res.has_value( j + 1 ) ) num_three_prime_connections++;

		bool const jump_to_next = check_jump_to_next_residue_in_chain( pose, j, slice_res );
		if ( jump_to_next ) num_jumps_to_next++;

		if ( rsd->is_polymer() && ( n == 1 || ( !after_cutpoint && slice_res.has_value( j - 1 ) ) ) ) {
			sliced_out_pose.append_residue_by_bond( *rsd, true /* build_ideal_geometry */ );
		} else {
			sliced_out_pose.append_residue_by_jump( *rsd, sliced_out_pose.size() );
		}
	}

	runtime_assert ( num_five_prime_connections <= 1 );
	runtime_assert ( num_three_prime_connections <= 1 );
	runtime_assert ( num_jumps_to_previous <= 1 );
	runtime_assert ( num_jumps_to_next <= 1 );
	// requirement for a clean slice:
	runtime_assert( (num_five_prime_connections + num_three_prime_connections + num_jumps_to_previous + num_jumps_to_next) == 1 );

	// fold tree!
	// figure out jumps
	utility::vector1< Size > jump_partners1, jump_partners2, endpoints, cuts;
	utility::vector1< std::string > jump_atoms1, jump_atoms2;
	FoldTree const & f = pose.fold_tree();
	for ( Size n = 1; n <= f.num_jump(); n++ ) {
		Size const j1 = f.upstream_jump_residue( n );
		Size const j2 = f.downstream_jump_residue( n );
		if ( slice_res.has_value( j1 ) && slice_res.has_value( j2 ) ) {
			jump_partners1.push_back( slice_res.index( j1 ) );
			jump_partners2.push_back( slice_res.index( j2 ) );
			jump_atoms1.push_back( f.upstream_atom( n ) );
			jump_atoms2.push_back( f.downstream_atom( n ) );
		}
	}

	for ( Size n = 1; n < slice_res.size(); n++ ) {
		Size const & i = slice_res[ n ];
		if ( !slice_res.has_value( i+1 ) || // boundary to non-sliced segment.
				pose.fold_tree().is_cutpoint( i ) ) {  // cut internal to sliced segment.
			cuts.push_back( n );
		}
	}

	Size const num_cuts = cuts.size();
	runtime_assert( jump_partners1.size() == jump_partners2.size() );
	runtime_assert( num_cuts == jump_partners1.size() );

	FoldTree f_slice;
	ObjexxFCL::FArray2D< int > jump_point_( 2, num_cuts, 0 );
	ObjexxFCL::FArray1D< int > cuts_( num_cuts, 0 );
	for ( Size i = 1; i <= num_cuts; i++ ) {
		jump_point_( 1, i ) = std::min( jump_partners1[ i ], jump_partners2[ i ] );
		jump_point_( 2, i ) = std::max( jump_partners1[ i ], jump_partners2[ i ] );
		cuts_( i ) = cuts[ i ];
	}
	f_slice.tree_from_jumps_and_cuts( sliced_out_pose.size(), num_cuts, jump_point_, cuts_ );

	// fix jump atoms.
	bool const KeepStubInResidue( true );
	for ( Size i = 1; i <= num_cuts; i++ ) {
		Size const n = f_slice.jump_nr( jump_partners1[ i ], jump_partners2[ i ] );
		f_slice.set_jump_atoms( n,
			jump_partners1[ i ], jump_atoms1[ i ],
			jump_partners2[ i ], jump_atoms2[ i ], KeepStubInResidue );
	}
	if ( slice_res.has_value( pose.fold_tree().root() ) ) f_slice.reorder( slice_res.index( pose.fold_tree().root() ) );

	sliced_out_pose.fold_tree( f_slice );

	// map (internal) coordinates from separate poses into merged on.
	// this is potentially dangerous, if we accumulate floating point errors, and scores change.
	std::map< Size, Size > res_map;
	for ( Size n = 1; n <= slice_res.size(); n++ ) res_map[ n ] =  slice_res[ n ];
	copy_dofs_match_atom_names( sliced_out_pose, pose, res_map, false /*backbone_only*/, false /*side_chain_only*/, false /*ignore_virtual*/ );
	align::superimpose_pose_legacy( sliced_out_pose, pose, res_map );

	declare_chemical_bonds_at_cutpoints( sliced_out_pose, pose, slice_res );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
check_jump_to_previous_residue_in_chain( pose::Pose const & pose, Size const i,
	utility::vector1< Size > const & current_element,
	FullModelInfo const & full_model_info ){
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & chains_in_full_model = get_chains_full( pose );
	return  check_jump_to_previous_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
check_jump_to_previous_residue_in_chain( pose::Pose const & pose, Size const i,
	utility::vector1< Size > const & current_element ){
	utility::vector1< Size > res_list, chains_in_full_model;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		res_list.push_back( n );
		chains_in_full_model.push_back( 1 ); // could easily fix this for virtual residue.
	}
	return  check_jump_to_previous_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
check_jump_to_previous_residue_in_chain( pose::Pose const & pose, Size const i,
	utility::vector1< Size > const & current_element,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & chains_in_full_model ){
	Size previous_res = i - 1;
	Size const current_chain = chains_in_full_model[ res_list[ i ] ];
	while ( previous_res >= 1 &&
			chains_in_full_model[ res_list[ previous_res ] ] == current_chain ) {
		if ( !current_element.has_value( res_list[ previous_res ] ) &&
				pose.fold_tree().jump_nr( previous_res, i ) > 0 ) return previous_res;
		previous_res--;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
check_jump_to_next_residue_in_chain( pose::Pose const & pose, Size const i,
	utility::vector1< Size > const & current_element,
	FullModelInfo const & full_model_info ){
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & chains_in_full_model = get_chains_full( pose );
	return  check_jump_to_next_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
check_jump_to_next_residue_in_chain( pose::Pose const & pose, Size const i,
	utility::vector1< Size > const & current_element ){
	utility::vector1< Size > res_list, chains_in_full_model;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		res_list.push_back( n );
		chains_in_full_model.push_back( 1 ); // could easily fix this for virtual residue.
	}
	return  check_jump_to_next_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
check_jump_to_next_residue_in_chain( pose::Pose const & pose, Size const i,
	utility::vector1< Size > const & current_element,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & chains_in_full_model ){
	Size subsequent_res = i + 1;
	Size const current_chain = chains_in_full_model[ res_list[ i ] ];
	while ( subsequent_res <= pose.size() &&
			chains_in_full_model[ res_list[ subsequent_res ] ] == current_chain ) {
		if ( !current_element.has_value( res_list[ subsequent_res ] ) &&
				pose.fold_tree().jump_nr( subsequent_res, i ) > 0 ) return subsequent_res;
		subsequent_res++;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
fix_up_jump_atoms( pose::Pose & pose ){
	// appears necessary for protein
	core::kinematics::FoldTree f = pose.fold_tree();
	f.reassign_atoms_for_intra_residue_stubs(); // it seems silly that we need to do this separately.
	pose.fold_tree( f );
}

/////////////////////////////////////////////////////////////////////////////////////////
void
fix_up_jump_atoms_and_residue_type_variants( pose::Pose & pose_to_fix ) {
	fix_up_jump_atoms( pose_to_fix );
	fix_up_residue_type_variants( pose_to_fix );
	update_constraint_set_from_full_model_info( pose_to_fix ); //atom numbers may have shifted around with variant changes.
	core::scoring::rna::clear_rna_scoring_info( pose_to_fix );
	protocols::rna::denovo::secstruct_legacy::clear_rna_secstruct_legacy_info( pose_to_fix );
}

/////////////////////////////////////////////////////////////////////////////////////////
// trying to track down a discrepancy in pose scores after switching pose focus...
void
check_scores_from_parts( std::string const & tag,
	Real const original_score,
	pose::Pose & pose,
	pose::PoseOP original_pose_clone,
	utility::vector1< PoseOP > const & other_pose_list,
	core::scoring::ScoreFunctionCOP scorefxn ) {

	using namespace core::scoring;

	TR << tag << std::endl;

	Real const pose_score = ( *scorefxn )( pose );
	scorefxn->show( pose );
	Real score = pose.energies().total_energies()[ intermol ] + pose.energies().total_energies()[ loop_close ];

	// "manual" check on each pose to track down discrepancies in scores after switching pose focus.
	ScoreFunctionOP scorefxn2 = scorefxn->clone();
	scorefxn2->set_weight( intermol, 0.0 );
	scorefxn2->set_weight( loop_close, 0.0 );
	TR << "EACH PART 0 " << std::endl;
	score += ( *scorefxn2 )( *original_pose_clone );
	scorefxn2->show( *original_pose_clone );

	for ( Size n = 1; n <= other_pose_list.size(); n++ ) {
		score += ( *scorefxn2 )( *( other_pose_list[ n ] ) );
		TR << "EACH PART " << n << std::endl;
		scorefxn2->show( *other_pose_list[ n ]  );
	}

	TR << TR.Magenta << tag << ":  ORIGINAL " << original_score << "  RESCORE " << pose_score << "  SUM OF PARTS " << score << TR.Reset << std::endl;

	core::scoring::loop_graph::LoopGraph loop_graph;
	loop_graph.update( pose, true /*verbose*/ );

}

/////////////////////////////////////////////////////////////////////////////////////////
void
switch_focus_to_other_pose( pose::Pose & pose,
	Size const & focus_pose_idx,
	core::scoring::ScoreFunctionCOP scorefxn /* = 0 */ )
{

	using namespace core::pose;
	using namespace core::pose::full_model_info;

	if ( focus_pose_idx == 0 ) return;

	Real const score_before_switch_focus = ( scorefxn != 0 ) ? (*scorefxn)( pose ) : 0.0;

	utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();

	// the original pose (or its clone, actually) is relieved of the
	// duty of holding information on other poses
	PoseOP original_pose_clone = pose.clone();
	utility::vector1< PoseOP > blank_pose_list;
	nonconst_full_model_info( *original_pose_clone ).set_other_pose_list( blank_pose_list );

	// debugging verbiage.
	// if ( scorefxn != 0 ) check_scores_from_parts(  "BEFORE_SWITCH_FOCUS", score_before_switch_focus, pose, original_pose_clone, other_pose_list, scorefxn );

	// need to shift focus to the other pose. It will now be responsible for holding
	// pointers to the other poses.
	PoseCOP other_pose = other_pose_list[ focus_pose_idx ];
	FullModelInfoOP new_full_model_info = const_full_model_info( *other_pose ).clone_info();

	// copy in pose list from original pose.
	utility::vector1< PoseOP > new_other_pose_list = new_full_model_info->other_pose_list();
	// for now -- but later could have pose 'trees', in which other_poses themselves have other_poses.
	runtime_assert( new_other_pose_list.size() == 0 );

	new_other_pose_list.push_back( original_pose_clone );
	for ( Size i = 1; i <= other_pose_list.size(); i++ ) {
		if ( i == focus_pose_idx ) continue;
		new_other_pose_list.push_back( other_pose_list[ i ]->clone() );
	}
	new_full_model_info->set_other_pose_list( new_other_pose_list );

	// OK, now shift focus! Hope this works.
	pose = ( *other_pose ); // makes a copy.
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, new_full_model_info );
	protocols::scoring::set_vdw_cached_rep_screen_info_from_pose( pose, *other_pose );

	Real const score_after_switch_focus = ( scorefxn != 0 ) ? (*scorefxn)( pose ) : 0.0;

	// debugging verbiage.
	// if ( scorefxn != 0 ) {
	//  PoseOP other_pose_clone = other_pose->clone();
	//  if ( scorefxn != 0 ) check_scores_from_parts(  "AFTER_SWITCH_FOCUS", score_after_switch_focus, pose, other_pose_clone, new_other_pose_list, scorefxn );
	// }

	if (  std::abs( score_before_switch_focus - score_after_switch_focus ) > 1.0e-3 ) {
		std::string const msg = "Energy change after switching pose focus: " + string_of( score_before_switch_focus ) + " to " +string_of( score_after_switch_focus );
		if ( !basic::options::option[ basic::options::OptionKeys::score::loop_close::allow_complex_loop_graph ]() ) {
			utility_exit_with_message( msg );
		}
		TR.Warning << TR.Red << msg << TR.Reset << std::endl;
	}


}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
switch_focus_among_poses_randomly( pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn /* = 0 */,
	bool force_switch /* = false */ ) {

	using namespace core::pose;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise;

	Size const num_other_poses = const_full_model_info( pose ).other_pose_list().size();
	if ( force_switch ) runtime_assert( num_other_poses > 0 );

	Size const focus_pose_idx = force_switch ?
		numeric::random::rg().random_range( 1, num_other_poses ) :
		numeric::random::rg().random_range( 0, num_other_poses );
	if ( focus_pose_idx == 0 ) return false;

	TR.Debug << TR.Green << "SWITCHING FOCUS! SWITCHING FOCUS! SWITCHING FOCUS! SWITCHING FOCUS! to: " << focus_pose_idx << TR.Reset << std::endl;
	switch_focus_to_other_pose( pose, focus_pose_idx, scorefxn );

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// specific for two helix test case. Sort of a unit test. Maybe I should make it a unit test.
void
test_merge_and_slice_with_two_helix_test_case(  utility::vector1< core::pose::PoseOP > const & input_poses,
	core::scoring::ScoreFunctionOP scorefxn ){

	using namespace core::pose;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise;

	scorefxn->show( *( input_poses[ 1 ] ) );
	scorefxn->show( *( input_poses[ 2 ] ) );

	input_poses[1]->dump_pdb( "input_pose1.pdb" );
	input_poses[2]->dump_pdb( "input_pose2.pdb" );

	pose::Pose merged_pose;
	merge_two_poses_using_full_model_info( merged_pose, *(input_poses[ 2 ]), *( input_poses[ 1 ] ),
		12 /*lower_merge_res*/, 13 /*upper_merge_res*/, true /*connect_by_bond*/ );
	merged_pose.dump_pdb( "merged_pose.pdb" );
	scorefxn->show( merged_pose );

	pose::Pose sliced_out_pose;
	utility::vector1< Size > slice_res = const_full_model_info( *(input_poses[2]) ).res_list();
	slice_out_pose( merged_pose, sliced_out_pose, slice_res );
	merged_pose.dump_pdb( "pose_after_slice.pdb" );
	sliced_out_pose.dump_pdb( "sliced_out_pose.pdb" );
	scorefxn->show( merged_pose );
	scorefxn->show( sliced_out_pose );
}


///////////////////////////////////////////////////////////////////
/// TODO: generalize to cyclized chainbreaks
utility::vector1< Size >
figure_out_moving_chain_break_res( pose::Pose const & pose, kinematics::MoveMap const & mm ) {

	utility::vector1< Size > moving_chainbreak_res;
	utility::vector1< Size > const cutpoint_closed = get_cutpoint_closed( pose );

	// kind of brute-force check over all moving residues. Can optimize if needed.
	for ( Size n = 1; n < pose.size(); n++ ) {

		if ( cutpoint_closed.has_value( n ) ) {
			// already a chainbreak -- is it moving?
			if ( moving_chainbreak_res.has_value( n ) ) continue;
			Size const n_next = core::scoring::methods::get_upper_cutpoint_partner_for_lower( pose, n );
			if ( mm.get( id::TorsionID( n, id::BB, 5 ) ) ||
					mm.get( id::TorsionID( n, id::BB, 6 ) ) ||
					mm.get( id::TorsionID( n_next, id::BB, 1 ) ) ||
					mm.get( id::TorsionID( n_next, id::BB, 2 ) ) ||
					mm.get( id::TorsionID( n_next, id::BB, 3 ) ) ) moving_chainbreak_res.push_back( n );
		} else {

			// not a chainbreak -- if we move it will any cutpoints move relative to each other?
			if ( pose.fold_tree().is_cutpoint( n ) ) continue;

			if ( ! mm.get( id::TorsionID( n+1, id::BB, 1 ) ) ) continue;

			ObjexxFCL::FArray1D_bool partner1( pose.size(), false );
			pose.fold_tree().partition_by_residue( n, partner1 );

			for ( Size m = 1; m <= cutpoint_closed.size(); m++ ) {
				Size const k = cutpoint_closed[ m ];
				if ( moving_chainbreak_res.has_value( k ) ) continue;
				Size const k_upper( core::scoring::methods::get_upper_cutpoint_partner_for_lower( pose, k ) );
				if ( partner1( k ) != partner1( k_upper ) ) moving_chainbreak_res.push_back( k );
			}
		}
	}

	return moving_chainbreak_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
check_for_input_domain( pose::Pose const & pose,
	utility::vector1< Size> const & partition_res ){
	utility::vector1< Size > const domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
	for ( Size i = 1; i <= partition_res.size(); i++ ) {
		if ( domain_map[ partition_res[ i ] ] > 0 ) return true;
	}
	return false;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
primary_fixed_domain( pose::Pose const & pose,
	utility::vector1< Size> const & partition_res ){
	utility::vector1< Size > const domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
	Size d_min( 0 );
	for ( Size i = 1; i <= partition_res.size(); i++ ) {
		Size const d = domain_map[ partition_res[ i ] ];
		if ( d > 0 && ( d_min == 0 || d < d_min ) ) d_min = d;
	}
	return d_min;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
check_for_input_domain( pose::Pose const & pose ){
	utility::vector1< Size > partition_res;
	for ( Size i = 1; i <= pose.size(); i++ ) partition_res.push_back( i );
	return check_for_input_domain( pose, partition_res );
}

///////////////////////////////////////////////////////////////
void
make_variants_match( pose::Pose & pose,
	pose::Pose const & reference_pose,
	Size const n,
	chemical::VariantType const variant_type )
{
	if ( reference_pose.residue_type( n ).has_variant_type( variant_type ) ) {
		add_variant_type_to_pose_residue( pose, variant_type, n );
	} else {
		remove_variant_type_from_pose_residue( pose, variant_type, n );
	}
}

//////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_moving_cutpoints_closed_from_moving_res( pose::Pose const & pose, Size const moving_res ){
	return figure_out_moving_cutpoints_closed( pose, figure_out_moving_partition_res( pose, moving_res ) );
}

//////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_moving_cutpoints_closed_from_moving_res( pose::Pose const & pose, utility::vector1< Size > const & moving_res_list ){
	utility::vector1< Size > moving_cutpoints;
	for ( Size n = 1; n <= moving_res_list.size(); n++ ) {
		moving_cutpoints = merge_vectors( moving_cutpoints,  figure_out_moving_cutpoints_closed( pose, figure_out_moving_partition_res( pose, moving_res_list[n] ) ) );
	}
	return moving_cutpoints;
}

//////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_moving_cutpoints_closed( pose::Pose const & pose,
	utility::vector1< Size > const & moving_partition_res ){
	utility::vector1< Size > cutpoints_closed, five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes;
	figure_out_moving_chain_breaks( pose, moving_partition_res, cutpoints_closed, five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes );
	return cutpoints_closed;
}

//////////////////////////////////////////////////////////////////////////////////
/// TODO: generalize to cyclized chainbreaks.
void
figure_out_moving_chain_breaks( pose::Pose const & pose,
	utility::vector1< Size > const & moving_partition_res,
	utility::vector1< Size > & cutpoints_closed,
	utility::vector1< Size > & five_prime_chain_breaks,
	utility::vector1< Size > & three_prime_chain_breaks,
	utility::vector1< Size > & chain_break_gap_sizes ){

	utility::vector1< Size > const chains = figure_out_chain_numbers_from_full_model_info_const( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	Size five_prime_chain_break( 0 ), three_prime_chain_break( 0 );
	for ( Size n = 1; n < pose.size(); n++ ) {

		if ( !pose.fold_tree().is_cutpoint( n ) ) continue;

		// Skip virtual anchors
		if ( pose.residue_type( n ).aa() == core::chemical::aa_vrt ) continue;
		if ( pose.residue_type( n+1 ).aa() == core::chemical::aa_vrt ) continue;

		// must be in different partitions to qualify as 'moving'
		if ( moving_partition_res.has_value( n ) == moving_partition_res.has_value( n+1 ) ) continue;

		// must be in same chain to qualify as a chain break
		if ( chains[ n ] != chains[ n+1 ] ) continue;

		// rewind to non-virtual residue closest to chainbreak
		for ( five_prime_chain_break = n; five_prime_chain_break >= 1; five_prime_chain_break-- ) {
			if ( !pose.residue_type( five_prime_chain_break ).has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) break;
		}

		// fast-forward to non-virtual residue closest to chainbreak
		for ( three_prime_chain_break = n+1; three_prime_chain_break <= pose.size(); three_prime_chain_break++ ) {
			if ( !pose.residue_type( three_prime_chain_break ).has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) break;
		}

		five_prime_chain_breaks.push_back( five_prime_chain_break );
		three_prime_chain_breaks.push_back( three_prime_chain_break );
		runtime_assert( res_list[ three_prime_chain_break ] > res_list[ five_prime_chain_break ] );
		Size const gap_size = res_list[ three_prime_chain_break ] - res_list[ five_prime_chain_break ] - 1;
		chain_break_gap_sizes.push_back( gap_size );
		if ( gap_size == 0 ) {
			runtime_assert( five_prime_chain_break == n ); // no rewind past bulges
			runtime_assert( three_prime_chain_break == n + 1 ); // no fast forward past bulges
			runtime_assert( pose.residue_type( five_prime_chain_break ).has_variant_type( chemical::CUTPOINT_LOWER ) );
			runtime_assert( pose.residue_type( three_prime_chain_break ).has_variant_type( chemical::CUTPOINT_UPPER ) );
			cutpoints_closed.push_back( n );
		}

	}

}

///////////////////////////////////////////////////////////////////////
utility::vector1< bool >
get_partition_definition_by_jump( pose::Pose const & pose, Size const & jump_nr /*jump_number*/ ){
	return pose.fold_tree().partition_by_jump( jump_nr );
}

///////////////////////////////////////////////////////////////////////
utility::vector1< bool >
get_partition_definition( pose::Pose const & pose, Size const & moving_suite ){

	ObjexxFCL::FArray1D<bool> partition_definition( pose.size(), false );
	if ( moving_suite > 0 ) {
		pose.fold_tree().partition_by_residue( moving_suite, partition_definition );
	} else {
		partition_definition.dimension( pose.size(), true );
	}

	//silly conversion. There may be a faster way to do this actually.
	utility::vector1< bool > partition_definition_vector1;
	for ( Size n = 1; n <= pose.size(); n++ ) partition_definition_vector1.push_back( partition_definition(n) );

	return partition_definition_vector1;

}


///////////////////////////////////////////////////////////////////////////
Size
figure_out_reference_res_for_suite( pose::Pose const & pose, Size const moving_res ){
	bool connected_by_jump( false );
	Size reference_res =  pose.fold_tree().get_parent_residue( moving_res, connected_by_jump );
	runtime_assert( !connected_by_jump );
	return reference_res;
}

///////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_moving_partition_res_for_suite( pose::Pose const & pose,
	Size const moving_res,
	Size const reference_res ){
	Size const moving_suite = ( moving_res < reference_res ) ? moving_res : reference_res;
	utility::vector1< bool > partition_definition = get_partition_definition( pose, moving_suite );
	utility::vector1< Size > moving_partition_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( partition_definition[ n ] == partition_definition[ moving_res ] ) moving_partition_res.push_back( n );
	}
	return moving_partition_res;
}

///////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_moving_partition_res_for_jump( pose::Pose const & pose,
	Size const jump_nr ){
	utility::vector1< bool > partition_definition = get_partition_definition_by_jump( pose, jump_nr );
	utility::vector1< Size > moving_partition_res;
	Size const moving_res = pose.fold_tree().downstream_jump_residue( jump_nr );
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( partition_definition[ n ] == partition_definition[ moving_res ] ) moving_partition_res.push_back( n );
	}
	return moving_partition_res;
}

///////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_moving_partition_res( pose::Pose const & pose,
	Size const moving_res ){
	utility::vector1< Size > root_partition_res, moving_partition_res;
	figure_out_root_and_moving_partition_res( pose, moving_res, root_partition_res, moving_partition_res );
	return moving_partition_res;
}

///////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_moving_partition_res( pose::Pose const & pose,
	utility::vector1< Size > const & moving_res_list ){

	utility::vector1< Size > moving_partition_res_all;
	for ( Size const moving_res : moving_res_list ) {
		if ( pose.residue_type( moving_res ).has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) ) continue; // legacy insanity from SWA 'dinucleotide' move.
		utility::vector1< Size > moving_partition_res = figure_out_moving_partition_res( pose, moving_res );
		for ( Size const seqpos : moving_partition_res ) {
			if ( !moving_partition_res_all.has_value( seqpos ) ) moving_partition_res_all.push_back( seqpos );
		}
	}
	std::sort( moving_partition_res_all.begin(), moving_partition_res_all.end() );
	return moving_partition_res_all;
}

///////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_root_partition_res( pose::Pose const & pose,
	utility::vector1< Size > const & moving_res_list ){
	utility::vector1< Size > root_partition_res;
	utility::vector1< Size > const moving_partition_res = figure_out_moving_partition_res( pose, moving_res_list );
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !moving_partition_res.has_value( n ) )  root_partition_res.push_back( n ) ;
	}
	return root_partition_res;
}

/////////////////////////////////////////////////////////////////////////////
void
figure_out_root_and_moving_partition_res( pose::Pose const & pose, Size const moving_res,
	utility::vector1< Size > & root_partition_res,
	utility::vector1< Size > & moving_partition_res ) {

	FoldTree const & f = pose.fold_tree();
	Size const reference_res = f.get_parent_residue( moving_res );
	Size const jump_nr = f.jump_nr( moving_res, reference_res );
	if ( !jump_nr && reference_res > 0 ) runtime_assert( moving_res == reference_res + 1 || reference_res == moving_res + 1 );
	utility::vector1< bool > partition_definition = jump_nr ? get_partition_definition_by_jump( pose, jump_nr ) : get_partition_definition( pose, std::min( moving_res, reference_res ) );
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( partition_definition[ n ] == partition_definition[ moving_res ] ) {
			moving_partition_res.push_back( n );
		} else {
			root_partition_res.push_back( n );
		}
	}
}

///////////////////////////////////////////////////////////////////
bool
revise_root_and_moving_res( pose::Pose & pose, Size & moving_res /* note that this can change too*/ ){
	if ( moving_res == 0 ) return false;

	utility::vector1< Size > root_partition_res, moving_partition_res;
	figure_out_root_and_moving_partition_res( pose, moving_res, root_partition_res, moving_partition_res );

	bool switch_moving_and_root_partitions( false );
	Size const primary_domain_moving = primary_fixed_domain( pose, moving_partition_res );
	Size const primary_domain_root   = primary_fixed_domain( pose, root_partition_res );
	if ( primary_domain_moving == 0 && primary_domain_root > 0 ) {
		switch_moving_and_root_partitions = false;
	} else if ( primary_domain_moving > 0 && primary_domain_root == 0 ) {
		switch_moving_and_root_partitions = true;
	} else if ( !get_definite_terminal_root( pose, moving_partition_res ) &&
			get_definite_terminal_root( pose, root_partition_res ) )  {
		switch_moving_and_root_partitions = false;
	} else if ( !get_definite_terminal_root( pose, root_partition_res ) &&
			get_definite_terminal_root( pose, moving_partition_res ) )  {
		switch_moving_and_root_partitions = true;
	} else {
		// either both are fixed or both are free.
		if ( root_partition_res.size() == moving_partition_res.size() ) {
			switch_moving_and_root_partitions = ( primary_domain_moving < primary_domain_root );
		} else {
			switch_moving_and_root_partitions = ( root_partition_res.size() < moving_partition_res.size() );
		}
	}

	if ( switch_moving_and_root_partitions ) { // the way things should be:
		Size const moving_res_original = moving_res; // to check switching went OK.
		Size const reference_res = pose.fold_tree().get_parent_residue( moving_res );
		moving_res = reference_res;
		reroot_based_on_full_model_info( pose, moving_partition_res /* new root_parition_res*/  );
		runtime_assert( pose.fold_tree().get_parent_residue( moving_res ) == static_cast<int>( moving_res_original ) );
	} else {
		reroot_based_on_full_model_info( pose, root_partition_res );
	}

	return switch_moving_and_root_partitions;
}


/////////////////////////////////////////////////////////////////////////////////////
bool
revise_root_and_moving_res_list( pose::Pose & pose,
	utility::vector1< Size > & moving_res_list /* note that this can change too*/ ){

	if ( moving_res_list.size() == 0 ) return false; // maybe after a delete -- just minimize.
	if ( pose.residue_type( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) return false; // for ERRASER.

	// weird case, happens in some protein early moves -- get the root out of the moving res!revise_root_and
	if ( moving_res_list.has_value( pose.fold_tree().root() ) ) {
		utility::vector1< Size > other_res = get_other_residues( moving_res_list, pose.size() );
		if ( other_res.size() == 0 ) return false; // from scratch.
		reroot_based_on_full_model_info( pose, other_res );
	}

	// find connection point to 'fixed res'
	Size moving_res_at_connection( 0 ), reference_res( 0 );
	for ( Size n = 1; n <= moving_res_list.size(); n++ ) {
		Size const & moving_res = moving_res_list[n];
		Size const parent_res = pose.fold_tree().get_parent_residue( moving_res );
		if ( moving_res_list.has_value( parent_res ) ) continue;
		runtime_assert( moving_res_at_connection == 0 ); // moving_res_list must be contiguous in fold tree.
		moving_res_at_connection = moving_res;
		reference_res            = parent_res;
		// break; // don't break, actually -- let's make sure there's only one connection point.
	}
	runtime_assert( moving_res_at_connection != 0 );

	if ( reference_res == 0 ) return false; // special case -- happens with 'from-scratch' moves.

	bool const is_jump = pose.fold_tree().jump_nr( reference_res, moving_res_at_connection ) > 0 ;
	bool switched_moving_and_root_partitions = revise_root_and_moving_res( pose, moving_res_at_connection );
	if ( is_jump && moving_res_list.size() != 1 ) runtime_assert( !switched_moving_and_root_partitions );

	// revise_root_and_moving_res() handled revision of single residue only... need to translate this
	// switch to the whole list of residues.
	if (  switched_moving_and_root_partitions ) {
		utility::vector1< Size > const moving_res_list_original = moving_res_list;
		moving_res_list = utility::tools::make_vector1( moving_res_at_connection );
		if ( is_jump ) { // easy switch-a-roo. jump connected single residues, and switch them.
			runtime_assert( moving_res_at_connection /*was switched*/  == reference_res );
		} else {
			// might be more residues. Example:
			//
			//              M   M   M   <-- moving_res_list (original)
			// ROOT...- 3 - 4 - 5 - 6 - 7 - ....
			//            M   M   M     <-- defines these moving connections
			//
			// after re-rooting, now should be:
			//
			//          M   M   M       <-- moving_res_list (new)
			//     ...- 3 - 4 - 5 - 6 - 7 - .... NEW ROOT
			//            M   M   M     <-- defines the same moving connections
			//
			for ( Size n = 1; n <= moving_res_list_original.size(); n++ ) {
				Size const & moving_res = moving_res_list_original[n];
				if ( moving_res_list_original.has_value( pose.fold_tree().get_parent_residue( moving_res ) ) ) moving_res_list.push_back( moving_res );
			}
		}
		runtime_assert( moving_res_list.size() == moving_res_list_original.size() );
	}

	return switched_moving_and_root_partitions;
}

/////////////////////////////////////////////////////////////////////////////////////
Size
split_pose( pose::Pose & pose, Size const moving_res, Size const reference_res ){

	Size jump_at_moving_suite = pose.fold_tree().jump_nr( moving_res, reference_res );
	if ( !jump_at_moving_suite ) {
		runtime_assert( (moving_res == reference_res + 1) || (moving_res == reference_res - 1) );
		Size const moving_suite = ( moving_res < reference_res ) ? moving_res : reference_res;
		runtime_assert( !pose.fold_tree().is_cutpoint( moving_suite ) );
		jump_at_moving_suite = make_cut_at_moving_suite( pose, moving_suite );
		if ( pose.residue_type( moving_suite + 1 ).is_RNA() ) add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, moving_suite+1 );
	}

	kinematics::Jump j = pose.jump( jump_at_moving_suite );
	j.set_translation( Vector( 1.0e4, 0.0, 0.0 ) );
	pose.set_jump( jump_at_moving_suite, j );
	return jump_at_moving_suite;
}

////////////////////////////////////////////////////////////////////////////////
void
split_pose( pose::Pose & pose, utility::vector1< Size > const & moving_res_list ){
	for ( Size n = 1; n <= moving_res_list.size(); n++ ) {
		Size const moving_res = moving_res_list[n];
		Size const reference_res = pose.fold_tree().get_parent_residue( moving_res );
		if ( reference_res == 0 ) continue; // happens when root atom is 'moving', e.g. in moves from_scratch.
		split_pose( pose, moving_res, reference_res );
	}
}

////////////////////////////////////////////////////////////////////////////////
void
fix_protein_jump_atom( pose::Pose & pose, Size const res, std::string const & atom_name ){
	using namespace core::kinematics;
	if ( !pose.residue_type( res ).is_protein() ) return;
	FoldTree f = pose.fold_tree();
	for ( Size n = 1; n <= f.num_jump(); n++ ) {
		if ( f.upstream_jump_residue( n ) == static_cast<int>( res ) ) {
			f.set_jump_atoms( n, atom_name, f.downstream_atom( n ), false /*intra_residue_stubs*/ );
		}
		if ( f.downstream_jump_residue( n ) == static_cast<int>( res ) ) {
			f.set_jump_atoms( n, f.upstream_atom( n ), atom_name, false /*intra_residue_stubs*/ );
		}
	}
	pose.fold_tree( f );
}


///////////////////////////////////////////////////////////////////////////////
void
add_to_pose_list( utility::vector1< core::pose::PoseOP > & pose_list, pose::Pose const & pose, std::string const & pose_tag ) {
	core::pose::PoseOP pose_op = pose.clone();
	tag_into_pose( *pose_op, pose_tag );
	pose_list.push_back( pose_op );
}

//////////////////////////////////////////////////////////////////////////////
// will deprecate this soon.
bool
is_protein( pose::Pose const & pose, utility::vector1< Size > const & moving_res_list ) {
	Size const example_res = ( moving_res_list.size() > 0 ) ? moving_res_list[1] : 1;
	if ( pose.residue_type( example_res ).is_protein() ) {
		return true;
	} else {
		runtime_assert( pose.residue_type( example_res ).is_RNA() );
		return false;
	}
}


///////////////////////////////////////////////////////////////////////////////////
// suites that connect different domains, e.g., different input structures.
utility::vector1< core::Size >
get_domain_boundary_suites( pose::Pose const & pose ){
	utility::vector1< Size > domain_boundary_suites;
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	for ( Size n = 1; n < pose.size(); n++ ) {
		if ( !cutpoint_open_in_full_model.has_value( res_list[ n ] ) &&
				!pose.fold_tree().is_cutpoint( n ) &&
				( res_list[ n + 1 ]  == res_list[ n ] + 1 ) &&
				( fixed_domain_map[ res_list[ n + 1 ] ] !=  fixed_domain_map[ res_list[ n ] ] ) &&
				!domain_boundary_suites.has_value( n ) ) {
			TR.Debug << "ADDING NEW SUITE TO BE MINIMIZED BASED ON LOCATION AT DOMAIN BOUNDARY: " << n << std::endl;
			domain_boundary_suites.push_back( n );
		}
	}
	return domain_boundary_suites;
}

/////////////////////////////////////////////////////////////////////////////////
// suites that connect different domains, e.g., different input structures. This
// function returns downstream residues within each suite connection.
utility::vector1< core::Size >
get_domain_boundary_res( pose::Pose const & pose ){

	// convert from domain boundaries suites to domain boundary residues.
	utility::vector1< core::Size > domain_boundary_suites = get_domain_boundary_suites( pose );

	utility::vector1< core::Size > domain_boundary_res;
	for ( Size n = 1; n <= domain_boundary_suites.size(); n++ ) {
		Size const & i = domain_boundary_suites[ n ];
		if ( pose.fold_tree().get_parent_residue( i ) == int(i+1) ) {
			domain_boundary_res.push_back( i );
		} else {
			runtime_assert( pose.fold_tree().get_parent_residue( i+1 ) == int(i) );
			domain_boundary_res.push_back( i+1 );
		}
	}
	return domain_boundary_res;
}

////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::Size >
get_moving_res_including_domain_boundaries( pose::Pose const & pose, utility::vector1< core::Size > const & moving_res_list ){
	return merge_vectors( moving_res_list, get_domain_boundary_res( pose ) );
}

////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_all_working_moving_res( working_parameters::StepWiseWorkingParametersCOP working_parameters) {
	return merge_vectors( working_parameters->working_moving_res_list(),
		working_parameters->working_bridge_res() );
}

////////////////////////////////////////////////////////////////////////////////
void
virtualize_side_chains( pose::Pose & pose ) {
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue_type( n ).is_RNA() ) {
			add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_O2PRIME_HYDROGEN, n );
		} else if ( pose.residue_type( n ).is_protein() ) {
			add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_SIDE_CHAIN, n );
		}
	}
}

//////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_all_residues( pose::Pose const & pose ){
	utility::vector1< Size > all_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		all_res.push_back( n );
	}
	return all_res;
}

//////////////////////////////////////////////////////////////////////
Size
find_downstream_connection_res( pose::Pose const & pose,
	utility::vector1< Size > const & moving_partition_res ){
	Size downstream_connection_res( 0 );
	for ( Size n = 1; n <= moving_partition_res.size(); n++ ) {
		if ( pose.residue_type( moving_partition_res[ n ] ).aa() == core::chemical::aa_h2o ) continue; // hack -- no waters
		Size const parent_res = pose.fold_tree().get_parent_residue( moving_partition_res[ n ] );
		if ( parent_res > 0 && !moving_partition_res.has_value( parent_res ) ) {
			if ( downstream_connection_res > 0 ) {
				TR << pose.annotated_sequence() << std::endl;
				TR << pose.fold_tree() << std::endl;
				TR << "downstream_connection_res already filled: " << downstream_connection_res << "  new downstream_connection_res " << moving_partition_res[n] << " goes with parent_res " << parent_res << std::endl;
			}
			runtime_assert( downstream_connection_res == 0 );
			downstream_connection_res = moving_partition_res[n];
		}
	}
	return downstream_connection_res;
}


//////////////////////////////////////////////////////////////////////
Size
get_unique_connection_res( pose::Pose const & pose,
	utility::vector1< Size > const & moving_partition_res ) {
	Size moving_res1 = find_downstream_connection_res( pose, moving_partition_res );

	utility::vector1< Size > const other_partition_res = get_other_residues( moving_partition_res, pose.size() );
	Size moving_res2 = find_downstream_connection_res( pose, other_partition_res );

	// can't have 'downstream connections' both ways -- something wrong with partition.
	runtime_assert( (moving_res1 == 0) || ( moving_res2 == 0 ) );
	if ( moving_res1 > 0 ) return moving_res1;
	if ( moving_res2 > 0 ) return moving_res2;
	return pose.fold_tree().root();
}


} //modeler
} //stepwise
} //protocols
