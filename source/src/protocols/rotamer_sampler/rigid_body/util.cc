// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/FloatingBaseUtil.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rotamer_sampler/rigid_body/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>
#include <numeric/xyz.functions.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.rigid_body.FloatingBaseUtil" );

using namespace core;
using namespace core::pose;

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_atom_coordinates( utility::vector1< std::pair < id::AtomID, numeric::xyzVector< Real > > > & xyz_list, Size const & seq_num,
												conformation::Residue const & rsd_at_origin, kinematics::Stub const & moving_res_base_stub ){
		xyz_list.clear();
		Vector atom_pos;
		for ( Size at = 1; at <= rsd_at_origin.natoms(); at++ ){
			id::AtomID const id( at, seq_num );
			get_specific_atom_coordinate( at, atom_pos, rsd_at_origin, moving_res_base_stub );
			xyz_list.push_back( std::make_pair( id, atom_pos ) );
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_specific_atom_coordinate( Size const atom_index,
																Vector & atom_pos,
																conformation::Residue const & rsd_at_origin,
																kinematics::Stub const & moving_res_base_stub ) {
		atom_pos = rsd_at_origin.xyz( atom_index );
		get_specific_atom_coordinate( atom_pos, moving_res_base_stub );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_specific_atom_coordinate( Vector & atom_pos,
																kinematics::Stub const & moving_res_base_stub ) {
		numeric::xyzVector< Real > const & new_centroid = moving_res_base_stub.v;
		numeric::xyzMatrix< Real > const & new_coordinate_matrix = moving_res_base_stub.M;
		atom_pos = new_coordinate_matrix * atom_pos; //I think the order here does matter.
		atom_pos = atom_pos + new_centroid; //I think the order here does matter.
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_specific_atom_coordinate( std::string const & atom_name,
															  Vector & atom_pos,
																conformation::Residue const & rsd_at_origin,
																kinematics::Stub const & moving_res_base_stub ) {
		get_specific_atom_coordinate( rsd_at_origin.atom_index( atom_name ), atom_pos, rsd_at_origin, moving_res_base_stub );
	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	figure_out_reference_res_for_jump( pose::Pose const & pose,
																		 Size const moving_res ){
		kinematics::FoldTree const & fold_tree = pose.fold_tree();
		for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
			if ( fold_tree.downstream_jump_residue( n ) == static_cast< int >( moving_res ) ) {
				Size const reference_res = fold_tree.upstream_jump_residue( n );
				return reference_res;
			}
		}
		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_partition_res( pose::Pose const & pose,
																	 Size const moving_res, Size const reference_res ){

		Size jump_number( 0 );
		kinematics::FoldTree const & fold_tree = pose.fold_tree();
		for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
			if ( fold_tree.downstream_jump_residue( n ) == static_cast< int >( moving_res ) &&
					fold_tree.upstream_jump_residue( n ) == static_cast< int >( reference_res ) ) {
				jump_number = n;
				break;
			}
		}

		ObjexxFCL::FArray1D<bool> partition_definition( pose.total_residue(), false );
		pose.fold_tree().partition_by_jump( jump_number, partition_definition );

		utility::vector1< Size > moving_partition_res;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( partition_definition( n ) == partition_definition( moving_res ) ) moving_partition_res.push_back( n );
		}
		return moving_partition_res;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Set sets the base to the origin....
	void
	set_to_origin( pose::Pose & pose, Size const seq_num,
								 utility::vector1< Size > const & move_res_list,
								 bool verbose = false ){

		using namespace chemical;
		using namespace scoring;
		using namespace pose;
		using namespace id;

		conformation::Residue const & rsd( pose.residue( seq_num ) );

		numeric::xyzVector< Real > centroid = chemical::rna::get_rna_base_centroid( rsd, verbose );
		numeric::xyzMatrix< Real > base_coordinate_matrix = chemical::rna::get_rna_base_coordinate_system( rsd, centroid );
		numeric::xyzMatrix< Real > invert_coordinate_matrix = inverse( base_coordinate_matrix );

		for ( Size n = 1; n <= move_res_list.size(); n++ ){
			Size const i = move_res_list[ n ];
			for ( Size at = 1; at <= pose.residue_type( i ).natoms(); at++ ){
				id::AtomID const id( at, i );
				pose.set_xyz( id, pose.xyz( id ) - centroid ); //I think the order here does matter. Translate centroid to origin.
				pose.set_xyz( id, invert_coordinate_matrix * pose.xyz( id ) ); //I think the order here does matter. Rotate coordinate so that it equal to Rosetta internal reference frame
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	set_to_origin( pose::Pose & pose, Size const seq_num, bool verbose ){
		utility::vector1< Size > const move_res_list = utility::tools::make_vector1( seq_num );
		set_to_origin( pose, seq_num,
									 move_res_list, verbose );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	kinematics::Stub
	initialize_stub( conformation::Residue const & rsd ){
		runtime_assert( rsd.is_RNA() ); // for now. can easily update in the future!
		scoring::rna::RNA_CentroidInfo rna_centroid_info;
		return rna_centroid_info.get_base_coordinate_system( rsd );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	PoseCOP
	transform_moving_partition_to_origin( pose::Pose const & pose_start,
																				Size const moving_res,
																				utility::vector1< Size > const & moving_partition_res ){
		pose::PoseOP pose = pose_start.clone();
		set_to_origin( *pose, moving_res, moving_partition_res );
		return pose;
	}


} //rigid_body
} //rotamer_sampler
} //protocols
