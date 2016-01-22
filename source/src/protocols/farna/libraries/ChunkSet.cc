// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/libraries/ChunkSet.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/libraries/ChunkSet.hh>
#include <protocols/farna/libraries/RNA_ChunkLibrary.hh> // for ROSETTA_LIBRARY_DOMAIN
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/toolbox/AtomID_Mapper.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/copydofs/util.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.libraries.ChunkSet" );

using namespace core::pose;

namespace protocols {
namespace farna {
namespace libraries {

///////////////////////////////////////////////////////////////////////
ChunkSet::ChunkSet(
	utility::vector1< core::pose::PoseOP > const & pose_list,
	ResMap const & res_map ):
	res_map_( res_map ),
	user_input_( true )
{
	// MiniPose is a more compact format than Pose.
	for ( Size n = 1; n <= pose_list.size(); n++ ) {
		mini_pose_list_.push_back( core::pose::MiniPoseOP( new core::pose::MiniPose( *(pose_list[n]) ) ) );
	}
	filter_poses_have_same_sequence_and_variants();
	setup_atom_id_mask_and_mapper( *( pose_list[ 1 ] ) );
}

///////////////////////////////////////////////////////////////////////
ChunkSet::ChunkSet(
	utility::vector1< core::pose::MiniPoseOP > const & mini_pose_list,
	core::pose::Pose const & example_pose,
	ResMap const & res_map ):
	mini_pose_list_( mini_pose_list ),
	res_map_( res_map ),
	user_input_( true )
{
	filter_poses_have_same_sequence_and_variants();
	setup_atom_id_mask_and_mapper( example_pose );
}

///////////////////////////////////////////////////////////////////////
ChunkSet::~ChunkSet() {}

///////////////////////////////////////////////////////////////////////
void
ChunkSet::setup_atom_id_mask_and_mapper( core::pose::Pose const & pose ) {
	setup_atom_id_mask( pose );
	setup_atom_id_mapper_to_vanilla_chunk_pose( pose );
}

///////////////////////////////////////////////////////////////////////
void
ChunkSet::setup_atom_id_mask( core::pose::Pose const & pose )
{
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {

		core::conformation::Residue rsd = pose.residue( i );
		for ( Size j = 1; j <= rsd.natoms(); j++ ) {
			atom_id_mask_[ core::id::AtomID( j, i ) ] = !rsd.is_virtual( j );
		}

		// special case for magnesium, which has a couple virtual atoms that need to get moved around and to define stubs.
		// not elegant, but I want to get this working.
		if ( rsd.name3() == " MG" ) {
			for ( Size j = 1; j <= rsd.natoms(); j++ )  atom_id_mask_[ core::id::AtomID( j, i ) ] = true;
		}
	}
}

///////////////////////////////////////////////////////////////////////
// @detailed vanilla means no variant_types. throughout RNA fragment assembly,
//  taking a convention that atom_id's should be mapped to such poses without variants
void
ChunkSet::setup_atom_id_mapper_to_vanilla_chunk_pose( core::pose::Pose const & pose ) {
	using namespace protocols::toolbox;
	atom_id_mapper_to_vanilla_chunk_pose_ = AtomID_MapperCOP( new AtomID_Mapper( pose, true /*map_to_vanilla_pose*/ ) );
}


///////////////////////////////////////////////////////////////////////
// upper_terminus & lower_terminus currently do not do anything to RNA.
utility::vector1< std::string >
remove_terminus_variant_types_for_rna( utility::vector1< std::string > const & types, char seq ) {
	utility::vector1< std::string > types_filtered;
	for ( Size n = 1; n <= types.size(); n++ ) {
		std::string const & type( types[ n ] );
		if ( core::chemical::rna::rna_nts.find( seq ) != std::string::npos &&
				( type == "UPPER_TERMINUS_VARIANT" || type == "LOWER_TERMINUS_VARIANT" ) ) continue;
		types_filtered.push_back( type );
	}
	return types_filtered;
}

///////////////////////////////////////////////////////////////////////
void
ChunkSet::filter_poses_have_same_sequence_and_variants()
{
	using namespace core::pose;

	utility::vector1< MiniPoseOP > filtered_mini_pose_list;
	runtime_assert( mini_pose_list_.size() > 0 );
	filtered_mini_pose_list.push_back( mini_pose_list_[ 1 ] );

	utility::vector1< std::string > fullname_list;
	std::vector< Size > oneletter_to_fullname_index;
	std::string one_letter_sequence;
	parse_sequence( mini_pose_list_[ 1 ]->sequence(), fullname_list, oneletter_to_fullname_index, one_letter_sequence );

	for ( Size n = 2; n <= mini_pose_list_.size(); n++ ) {
		MiniPose const & mini_pose1 = *mini_pose_list_[ 1 ];
		MiniPose const & mini_pose2 = *mini_pose_list_[ n ];
		runtime_assert( mini_pose1.sequence() == mini_pose2.sequence() );
		for ( Size m = 1; m <= mini_pose1.total_residue(); m++ ) {
			utility::vector1< std::string > const & types1 = remove_terminus_variant_types_for_rna( mini_pose1.variant_types( m ), one_letter_sequence[m-1] );
			utility::vector1< std::string > const & types2 = remove_terminus_variant_types_for_rna( mini_pose2.variant_types( m ), one_letter_sequence[m-1] );
			if ( types1 != types2 ) {
				static Size count( 0 );
				count++;
				if ( count < 5 ) {
					TR << "filtering out of chunk_set pose with sequence " << mini_pose1.sequence() << " due to mismatch in variants at position " << m << ":  " <<  types1 << " vs " << types2 << TR.Reset << std::endl;
				} else if ( count == 5 ) {
					TR << "More chunk_set poses filtered out due to variant mismatch... suppressing these warnings." << std::endl;
				}
				continue;
			}
			filtered_mini_pose_list.push_back( mini_pose_list_[ n ] );
		}
	}
	mini_pose_list_ = filtered_mini_pose_list;
}


///////////////////////////////////////////////////////////////////////
void
ChunkSet::insert_chunk_into_pose( core::pose::Pose & pose, Size const & chunk_pose_index,
	toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	bool do_rosetta_library_domain_check /* = true */ ) const
{

	using namespace core::pose;
	using namespace core::id;

	core::pose::MiniPose const & scratch_pose ( *(mini_pose_list_[ chunk_pose_index ]) );

	std::map< AtomID, AtomID > atom_id_map = get_atom_id_map( pose, *atom_level_domain_map->atom_id_mapper() );
	std::map< AtomID, Size > atom_id_domain_map;
	if ( !user_input() ) atom_id_domain_map = get_atom_id_domain_map_for_rosetta_library_chunk( atom_id_map, pose, *atom_level_domain_map, do_rosetta_library_domain_check );

	core::pose::copydofs::copy_dofs( pose, scratch_pose,
		atom_id_map, atom_id_domain_map );

}

//////////////////////////////////////////////////////////////////////////////////////////////
std::map< id::AtomID, id::AtomID >
ChunkSet::get_atom_id_map(  core::pose::Pose & pose, toolbox::AtomID_Mapper const & atom_id_mapper_to_target_vanilla_pose ) const{

	std::map< id::AtomID, id::AtomID > atom_id_map = atom_id_mapper_to_target_vanilla_pose.calculate_atom_id_map( pose, res_map_,
		mini_pose_list_[1]->fold_tree(),
		atom_id_mapper_to_vanilla_chunk_pose_  );

	// This should prevent copying dofs for virtual phosphates, if they are tagged as such in the input silent files.
	filter_atom_id_map_with_mask( atom_id_map );

	return atom_id_map;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
ChunkSet::filter_atom_id_map_with_mask( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ) const{

	using namespace core::id;

	std::map< AtomID, AtomID > atom_id_map_new;

	for ( std::map< AtomID, AtomID >::const_iterator
			it=atom_id_map.begin(), it_end = atom_id_map.end(); it != it_end; ++it ) {

		AtomID const & target_atom_id = it->first;
		AtomID const & source_atom_id = it->second;

		std::map< AtomID, bool >::const_iterator it_mask = atom_id_mask_.find( source_atom_id );

		if ( it_mask == atom_id_mask_.end() ) utility_exit_with_message( "Some problem with atom_id_mask in defining atom_id_map " );

		if ( !it_mask->second ) continue;

		atom_id_map_new[ target_atom_id ] = source_atom_id;
	}

	atom_id_map = atom_id_map_new;

}

//////////////////////////////////////////////////////////////////////////////////////////////
// 'Rosetta library' domains are not user input, and they should not be put into
// the pose if some atoms are user inputted. At this point, atom_id_map should only
// contain atoms whose dofs might be affected by copy_dofs (filter_atom_id_with_mask better have
// been run). So the atoms either contain user-inputted domains or rosetta domains.
// This function sets up an atom-wise domain map with 0 at rosetta domains (moveable)
//////////////////////////////////////////////////////////////////////////////////////////////
std::map< core::id::AtomID, core::Size >
ChunkSet::get_atom_id_domain_map_for_rosetta_library_chunk(
	std::map< id::AtomID, id::AtomID > atom_id_map,
	pose::Pose const & pose, toolbox::AtomLevelDomainMap const & atom_level_domain_map,
	bool do_rosetta_library_domain_check /* = true */ ) const
{
	using namespace core::id;
	runtime_assert( !user_input() ); // we are a in rosetta library ChunkSet, not a user-inputted ChunkSet.

	std::map< AtomID, Size > atom_id_domain_map;
	bool found_rosetta_library_domain( false );

	for ( std::map< AtomID, AtomID >::const_iterator it=atom_id_map.begin(),
			it_end = atom_id_map.end(); it != it_end; ++it ) {
		AtomID const & target_atom_id = it->first;
		Size domain( atom_level_domain_map.get_domain( target_atom_id ) );
		if ( domain == ROSETTA_LIBRARY_DOMAIN ) {
			atom_id_domain_map[ target_atom_id ] = 0; // OK to insert here.
			found_rosetta_library_domain = true;
		} else {
			if ( do_rosetta_library_domain_check && domain == 0 ) {
				/// following is verbiage helpful for debugging. Remove after 2015 if not in use.
				atom_level_domain_map.show();
				for ( std::map< AtomID, AtomID >::const_iterator itx=atom_id_map.begin(),
						itx_end = atom_id_map.end(); itx != itx_end; ++itx ) {
					TR << TR.Green << atom_id_to_named_atom_id( itx->first, pose ) << " in pose with domain " << atom_level_domain_map.get_domain( itx->first ) << " mapped to " << itx->second << " in chunk " << std::endl;
				}

				TR << "mini pose has sequence: " << mini_pose_list_[1]->sequence() << std::endl;
				for ( Size m = 1; m <= mini_pose_list_[1]->total_residue(); m++ ) {
					TR << "In chunk mini-pose, residue " << m << " has variants: ";
					TR << mini_pose_list_[1]->variant_types( m ) << std::endl;
				}
				TR << "mini pose has fold tree: " << mini_pose_list_[1]->fold_tree() << std::endl;
				TR << pose.annotated_sequence() << std::endl;
				/// end verbiage.
				utility_exit_with_message( "Atom " + atom_id_to_named_atom_id( target_atom_id, pose ).to_string() + " should have domain > 0." );
			}
			atom_id_domain_map[ target_atom_id ] = domain; // not OK to insert here.
		}
	}
	if ( do_rosetta_library_domain_check ) runtime_assert( found_rosetta_library_domain );
	return atom_id_domain_map;
}


//////////////////////////////////////////////////////////////////////////////////////////////
core::pose::MiniPoseOP const
ChunkSet::mini_pose( Size const idx ) const {
	return mini_pose_list_[ idx ];
}

//////////////////////////////////////////////////////////////////////////////
bool
ChunkSet::check_fold_tree_OK( pose::Pose const & pose ) const
{
	// Check where the chunk is mapped to in the big pose.
	// There should be at least the same number of jumps in the big pose
	//  as there are chains in the scratch_pose.
	utility::vector1< bool > is_chunk_res( pose.total_residue(), false );
	for ( ResMap::const_iterator
			it=res_map_.begin(), it_end = res_map_.end(); it != it_end; ++it ) {
		Size const i = it->first; //Index in big pose.
		is_chunk_res[ i ] = true;
	}

	Size const num_jumps_scratch = mini_pose_list_[1]->fold_tree().num_jump(); // number of chains - 1

	Size num_jumps_in_big_pose_in_scratch_region( 0 );
	for ( Size n = 1; n <= pose.num_jump(); n++ ) {
		if ( ! is_chunk_res[ pose.fold_tree().upstream_jump_residue( n ) ] ) continue;
		if ( ! is_chunk_res[ pose.fold_tree().downstream_jump_residue( n ) ] ) continue;
		num_jumps_in_big_pose_in_scratch_region++;
	}

	if ( num_jumps_scratch > num_jumps_in_big_pose_in_scratch_region ) {
		std::cout << "Number of jumps in chunk pose               : " << num_jumps_scratch << std::endl;
		std::cout << "Number of jumps in full pose in chunk region: " << num_jumps_in_big_pose_in_scratch_region  << "  out of total jumps " << pose.num_jump() << std::endl;
		return false;
	}

	if ( num_jumps_scratch < num_jumps_in_big_pose_in_scratch_region ) {
		//   std::cout << "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!" << std::endl;
		//   std::cout << "Number of jumps in chunk pose               : " << num_jumps_scratch << std::endl;
		//   std::cout << "Does not match:" << std::endl;
		//   std::cout << "Number of jumps in full pose in chunk region: " << num_jumps_in_big_pose_in_scratch_region  << "  out of total jumps " << pose.num_jump() << std::endl;
		//   std::cout << "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!" << std::endl;
		// Just a warning
		//return false;
	}

	return true;

}

} //libraries
} //farna
} //protocols
