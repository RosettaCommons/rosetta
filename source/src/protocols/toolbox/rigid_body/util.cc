// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/rigid_body/FloatingBaseUtil.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/toolbox/rigid_body/util.hh>
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
#include <numeric/conversions.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.rigid_body.FloatingBaseUtil" );

using namespace core;
using namespace core::pose;
using numeric::conversions::degrees;
using numeric::conversions::radians;

using utility::operator <<;

//////////////////////////////////////////////////////////////////
//
//  Euler rotation utils developed in 2011-2012 for enumerative
//    rigid body sampling for, e.g., RNA nearest-neighbor rules.
//
//  Now mostly deprecated in favor of
//      stepwise/sampler/RigidBodyStepWiseSampler.hh
//
//  Just in use in nucleobase_sample_around.
//
//                               -- rhiju, 2014
//
//////////////////////////////////////////////////////////////////

namespace protocols {
namespace toolbox {
namespace rigid_body {

void
get_atom_coordinates( utility::vector1< std::pair < id::AtomID, numeric::xyzVector< Real > > > & xyz_list, Size const & seq_num,
	conformation::Residue const & rsd_at_origin, kinematics::Stub const & moving_res_base_stub ){
	xyz_list.clear();
	Vector atom_pos;
	for ( Size at = 1; at <= rsd_at_origin.natoms(); at++ ) {
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
		if ( fold_tree.downstream_jump_residue( n ) == moving_res ) {
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
		if ( fold_tree.downstream_jump_residue( n ) == moving_res &&
				fold_tree.upstream_jump_residue( n ) == reference_res ) {
			jump_number = n;
			break;
		}
	}

	ObjexxFCL::FArray1D<bool> partition_definition( pose.size(), false );
	pose.fold_tree().partition_by_jump( jump_number, partition_definition );

	utility::vector1< Size > moving_partition_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
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

	for ( Size n = 1; n <= move_res_list.size(); n++ ) {
		Size const i = move_res_list[ n ];
		for ( Size at = 1; at <= pose.residue_type( i ).natoms(); at++ ) {
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

	runtime_assert( rsd.is_RNA() || rsd.type().name() == "pdb_GAI" ); // for now. can easily update in the future!
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

///////////////////////////////////////////////////////////////////////
void
create_euler_rotation(
	Matrix & M,
	Real const & alpha,
	Real const & beta,
	Real const & gamma,
	Vector const & /* axis1 not actually used*/,
	Vector const & axis2,
	Vector const & axis3
)
{
	// Z-axis assumed to be long axis.
	Matrix M1 = numeric::rotation_matrix( axis3, Real( radians( alpha ) ) );
	Matrix M2 = numeric::rotation_matrix( axis2, Real( radians( beta ) ) );
	Matrix M3 = numeric::rotation_matrix( axis3, Real( radians( gamma ) ) );

	M = M3 * M2 * M1;
}

///////////////////////////////////////////////////////////////////////
void
create_euler_rotation(
	Matrix & M,
	Real const & alpha,
	Real const & beta,
	Real const & gamma )
{
	static Vector const axis1( 1.0, 0.0, 0.0 );
	static Vector const axis2( 0.0, 1.0, 0.0 );
	static Vector const axis3( 0.0, 0.0, 1.0 );
	create_euler_rotation( M, alpha, beta, gamma, axis1, axis2, axis3 );
}

//////////////////////////////////////////////////////////////////////////////////
void
get_euler_angles( Real & alpha, Real & beta, Real & gamma, Matrix M1, Matrix M2, bool const verbose /*=true*/ ){

	Matrix M_test = M2;

	// Figure out what axis system2 looks like in axis system1.
	M2 = M1.transposed() * M2;

	// First figure out how to backrotate z rotation.
	Vector z_vec = M2.col_z();
	Real const gamma_radians = std::atan2( z_vec(2), z_vec(1) );
	gamma = degrees( gamma_radians );

	M2 = rotation_matrix( Vector( 0.0, 0.0, 1.0 ), -1.0 * gamma_radians ) * M2;
	z_vec = M2.col_z();
	if ( verbose ) std::cout << "This better have a zero in y position " << z_vec(1) << ' ' << z_vec(2) << ' ' << z_vec(3) << std::endl;

	// Then figure out how to backrotate y rotation.
	Real const beta_radians = std::atan2( z_vec(1), z_vec(3) );
	beta = degrees( beta_radians );

	M2 = rotation_matrix( Vector( 0.0, 1.0, 0.0 ), -1.0 * beta_radians ) * M2;
	z_vec = M2.col_z();
	if ( verbose ) std::cout << "This better have a zero in x and y position " << z_vec(1) << ' ' << z_vec(2) << ' ' << z_vec(3) << std::endl;

	// Finally, backrotate z rotation.
	Vector x_vec = M2.col_x();
	Real const alpha_radians = std::atan2( x_vec(2), x_vec(1) );
	alpha = degrees( alpha_radians );

	M2 = rotation_matrix( Vector( 0.0, 0.0, 1.0 ), -1.0 * alpha_radians ) * M2;
	x_vec = M2.col_x();
	if ( verbose ) std::cout << "This better have a zero in y and z position " << x_vec(1) << ' ' << x_vec(2) << ' ' << x_vec(3) << std::endl;

	if ( verbose ) {
		Matrix M;

		Vector const xaxis1 = M1.col_x();
		Vector const yaxis1 = M1.col_y();
		Vector const zaxis1 = M1.col_z();
		create_euler_rotation( M, alpha, beta, gamma, xaxis1, yaxis1, zaxis1 );

		// Matrix M_test;
		// create_euler_rotation( M_test, alpha, beta, gamma, Vector( 1.0,0.0,0.0), Vector( 0.0,1.0,0.0), Vector(0.0,0.0,1.0) );
		// M_test = M1 * M_test * M1.transposed();
		// M_test = M1.transposed() * M_test * M1; // Can we rotate M1 into M2?
		M = M * M1;

		std::cout << "These better match:  " << std::endl;
		std::cout << M(1,1) <<  ' ' << M(1,2)  << ' ' << M(1,3) << std::endl;
		std::cout << M(2,1) <<  ' ' << M(2,2)  << ' ' << M(2,3) << std::endl;
		std::cout << M(3,1) <<  ' ' << M(3,2)  << ' ' << M(3,3) << std::endl;
		std::cout << std::endl;
		std::cout << M_test(1,1) <<  ' ' << M_test(1,2) <<  ' ' << M_test(1,3) << std::endl;
		std::cout << M_test(2,1) <<  ' ' << M_test(2,2)  << ' ' << M_test(2,3) << std::endl;
		std::cout << M_test(3,1) <<  ' ' << M_test(3,2)  << ' ' << M_test(3,3) << std::endl;
	}

}

///////////////////////////////////////////////////////////////////////
void
translate( pose::Pose & pose, Vector const & shift,
	pose::Pose const & ref_pose,
	utility::vector1< Size > const & moving_res ){

	using namespace core::id;

	for ( Size n = 1; n <= moving_res.size(); n++ ) {
		Size const i = moving_res[ n ];

		for ( Size m = 1; m <= pose.residue_type( i ).natoms(); m++ ) {
			pose.set_xyz( AtomID(m,i),   ref_pose.xyz( AtomID(m,i) ) + shift );
		}

	}

}


///////////////////////////////////////////////////////////////////////
void
rotate( pose::Pose & pose, Matrix const & M,
	pose::Pose const & ref_pose,
	utility::vector1< Size > const & moving_res,
	Vector const & centroid ){

	using namespace core::id;

	for ( Size n = 1; n <= moving_res.size(); n++ ) {
		Size const i = moving_res[ n ];

		for ( Size m = 1; m <= pose.residue_type( i ).natoms(); m++ ) {
			pose.set_xyz( AtomID(m,i),  M * ( ref_pose.xyz( AtomID(m,i) ) - centroid ) + centroid );
		}
	}

}

//////////////////////////////////////////////////////////////////
void
rotate( pose::Pose & pose, Matrix const & M,
	pose::Pose const & ref_pose,
	utility::vector1< Size > const & moving_res ){

	Vector centroid( 0.0, 0.0, 0.0 );
	rotate( pose, M, ref_pose, moving_res, centroid );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
get_base_centroid_and_rotation_matrix( pose::Pose const & pose, Size const i, Vector & centroid, Matrix & M ){

	using namespace scoring::rna;
	using namespace kinematics;
	static RNA_CentroidInfo rna_centroid_info;

	centroid = rna_centroid_info.get_base_centroid( pose.residue( i ) );
	Stub s = rna_centroid_info.get_base_coordinate_system( pose.residue( i ), centroid );
	M = s.M;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i, utility::vector1< Size > const & moving_res, bool const do_not_rotate  ){

	Vector centroid;
	Matrix M;
	get_base_centroid_and_rotation_matrix( pose, i, centroid, M);

	translate( pose, -centroid, pose, moving_res);
	if ( !do_not_rotate ) rotate( pose, M.transposed(), pose, moving_res);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i ){
	translate_and_rotate_residue_to_origin( pose, i, utility::tools::make_vector1( i ) );
}


} //rigid_body
} //toolbox
} //protocols
