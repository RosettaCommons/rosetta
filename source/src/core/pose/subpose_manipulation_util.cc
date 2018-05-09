// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/subpose_manipulation_util.cc
/// @brief  Pose class utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Rhiju Das, Steven Lewis, Vikram K. Mulligan


// Unit header
#include <core/pose/subpose_manipulation_util.hh>

// Package headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>

// Project headers
#include <core/chemical/rna/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

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

static basic::Tracer TR( "core.pose.subpose_manipulation_util" );

void
append_pose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	bool new_chain
){
	append_subpose_to_pose(pose1, pose2, 1, pose2.size(), new_chain);
}

void
append_subpose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::Size start_res,
	core::Size end_res,
	bool new_chain
){
	if ( pose2.size()<start_res ) {
		TR.Error << "Provided starting residue number " << start_res
			<< " less than number residues in appended pose. Nothing to do." << std::endl;
	}
	pose1.append_residue_by_jump(pose2.residue(start_res), pose1.size() , "", "", new_chain);
	for ( core::Size i=start_res+1; i<=end_res; ++i ) {
		if ( pose2.residue(i).is_lower_terminus() ) {
			if ( i > 1 && pose2.chain(i) == pose2.chain(i-1) ) {
				pose1.append_residue_by_jump(pose2.residue(i), pose1.size(), "","", false);
			} else {
				if ( pose2.residue(i).is_protein() ) {
					pose1.append_residue_by_bond(pose2.residue(i));
				} else {
					bool new_chain =  pose2.chain(i) != pose2.chain(i-1);
					pose1.append_residue_by_jump(pose2.residue(i), i-1, "","", new_chain);
				}
			}
		} else {
			pose1.append_residue_by_bond(pose2.residue(i));
		}
	}
}

void
create_subpose(
	Pose const & src,
	utility::vector1< Size > const & positions,
	kinematics::FoldTree const & f,
	Pose & pose
) {
	Size const nres( f.nres() );
	debug_assert( nres == positions.size() );

	pose.clear();

	for ( Size i=1; i<= nres; ++i ) {
		Size const seqpos( positions[i] );
		conformation::Residue const & rsd( src.residue( seqpos ) );
		// If the residue and the previous residue are bonded in the source pose, they should be bonded in the new pose
		if ( i>1 && rsd.is_polymer_bonded( positions[ i-1 ] ) ) {
			pose.append_residue_by_bond( rsd );
		} else {
			pose.append_residue_by_jump( rsd, 1 );
		}
		if ( i>1 ) {
			// check if this residue should be in a new chain. not a perfect check...
			conformation::Residue const & prev_rsd( src.residue( positions[i-1] ) );
			if ( prev_rsd.is_upper_terminus() || rsd.is_lower_terminus() || prev_rsd.chain() != rsd.chain() ) {
				debug_assert( pose.size() == i );
				pose.conformation().insert_chain_ending( i-1 );
			}
		}
	}

	// now set the desired foldtree
	pose.fold_tree(f);
}


///////////////////////////////////////////////////////////////////////////////
// A bit like create_subpose() but does not require fold tree, and
// a little optimized to handle some weird terminal variants and RNA jumps.
void
pdbslice( core::pose::Pose & new_pose,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & slice_res )
{
	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::chemical::rna;

	new_pose.clear();

	for ( Size i = 1; i <= slice_res.size(); i++ ) {

		ResidueOP residue_to_add = pose.residue( slice_res[ i ] ).clone() ;

		if ( (i > 1 &&  ( slice_res[i] != slice_res[i-1] + 1 )) /*new segment*/ ||
				residue_to_add->is_lower_terminus() ||
				residue_to_add->has_variant_type( N_ACETYLATION ) ||
				!residue_to_add->is_polymer() ||
				(i>1 && pose.fold_tree().is_cutpoint( slice_res[i-1] ) ) ) {

			if ( residue_to_add->is_RNA() && (i>1) && new_pose.residue_type(i-1).is_RNA() ) {

				new_pose.append_residue_by_jump(  *residue_to_add, i-1,
					chi1_torsion_atom( new_pose.residue_type(i-1) ),
					chi1_torsion_atom( residue_to_add->type() ), true /*new chain*/ );
			} else {

				new_pose.append_residue_by_jump(  *residue_to_add, i-1, "", "", true /*new chain*/ );
			}
		} else {

			new_pose.append_residue_by_bond(  *residue_to_add  ) ;
		}
	}

	using namespace core::pose::full_model_info;
	if ( full_model_info_defined( pose ) ) {
		FullModelInfoOP full_model_info = const_full_model_info( pose ).clone_info();
		utility::vector1< Size > const & res_list = full_model_info->res_list();
		utility::vector1< Size > new_res_list;
		for ( Size n = 1; n <= new_pose.size(); n++ ) {
			new_res_list.push_back( res_list[ slice_res[ n ] ] );
		}
		full_model_info->set_res_list( new_res_list );
		set_full_model_info( new_pose, full_model_info );
	}

	PDBInfoCOP pdb_info = pose.pdb_info();
	if ( pdb_info ) {
		utility::vector1< Size > new_numbering;
		utility::vector1< char > new_chains;
		utility::vector1< std::string > new_segids;
		for ( Size n = 1; n <= slice_res.size(); n++ ) {
			new_numbering.push_back( pdb_info->number( slice_res[ n ] ) );
			new_chains.push_back( pdb_info->chain( slice_res[ n ] ) );
			new_segids.push_back( pdb_info->segmentID( slice_res[ n ] ) );
		}
		PDBInfoOP new_pdb_info( new PDBInfo( new_pose, true /*init*/ ) );
		new_pdb_info->set_numbering( new_numbering );
		new_pdb_info->set_chains( new_chains );
		new_pdb_info->set_segment_ids( new_segids );
		new_pose.pdb_info( new_pdb_info );
	}
	tag_into_pose( new_pose, tag_from_pose( pose ) );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pdbslice(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & slice_res
) {
	// If slice_res is 'every residue in pose' just return
	// This happens in rna_denovo where working_res is empty (then
	// filled to 'everything' in creating the working_native.
	utility::vector1< Size > identity_slice_res( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) identity_slice_res[ ii ] = ii;
	if ( slice_res == identity_slice_res ) return;

	pose::Pose mini_pose;
	pdbslice( mini_pose, pose, slice_res );
	pose = mini_pose;
}

////////////////////////////////////////////////////////////////////////////
void
partition_pose_by_jump(
	pose::Pose const & src,
	int const jump_number,
	pose::Pose & partner1,
	pose::Pose & partner2
) {
	Size const nres( src.size() );

	// split src pose's foldtree
	kinematics::FoldTree f1, f2;
	src.fold_tree().partition_by_jump( jump_number, f1, f2 );

	TR << src.fold_tree() << '\n' << f1 << '\n' << f2 << std::endl;

	// identify residues in the two partners
	ObjexxFCL::FArray1D_bool partner1_pos( nres, false ); // FARRAY! DOH!!
	src.fold_tree().partition_by_jump( jump_number, partner1_pos );

	utility::vector1< Size > partner1_pos_list, partner2_pos_list;
	for ( Size i=1; i<= nres; ++i ) {
		if ( partner1_pos(i) ) partner1_pos_list.push_back( i );
		else partner2_pos_list.push_back( i );
	}

	create_subpose( src, partner1_pos_list, f1, partner1 );
	create_subpose( src, partner2_pos_list, f2, partner2 );
}


} // pose
} // core
