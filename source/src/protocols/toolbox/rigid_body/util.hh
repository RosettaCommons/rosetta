// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/rigid_body/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_toolbox_rigid_body_FloatingBaseUtil_HH
#define INCLUDED_protocols_toolbox_rigid_body_FloatingBaseUtil_HH

#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzMatrix.hh>
#include <string>

using namespace core;
using namespace core::pose;

typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace toolbox {
namespace rigid_body {

void
get_atom_coordinates( utility::vector1< std::pair < id::AtomID, numeric::xyzVector< Real > > > & xyz_list,
	Size const & seq_num,
	conformation::Residue const & rsd_at_origin,
	kinematics::Stub const & moving_res_base_stub );

void
get_specific_atom_coordinate( Vector & atom_pos,
	kinematics::Stub const & moving_res_base_stub );

void
get_specific_atom_coordinate( Size const atom_index,
	Vector & atom_pos,
	conformation::Residue const & rsd_at_origin,
	kinematics::Stub const & moving_res_base_stub );


void
get_specific_atom_coordinate( std::string const & atom_name,
	Vector & atom_pos,
	conformation::Residue const & rsd_at_origin,
	kinematics::Stub const & moving_res_base_stub );


/////////////////////////////////////////////////////////////
// following functions are helpers for RigidBodyStepWiseSampler
Size
figure_out_reference_res_for_jump( pose::Pose const & pose,
	Size const moving_res );

utility::vector1< Size >
figure_out_moving_partition_res( pose::Pose const & pose,
	Size const moving_res, Size const reference_res );

void
set_to_origin( pose::Pose & pose, Size const seq_num, bool verbose );

kinematics::Stub
initialize_stub( conformation::Residue const & rsd );

PoseCOP
transform_moving_partition_to_origin( pose::Pose const & pose_start,
	Size const moving_res,
	utility::vector1< Size > const & moving_partition_res_ );

/////////////////////////////////////////////////////////////////////////////////////////
// Following were developed for a rigid-body rotation enumeration framework for
// computing partition functions for the RNA nearest-neighbor rules. No
// longer in use. However, could be very useful in lo-res sampling, so I didn't
// remove. If these are not in use by 2016, delete!
void
get_euler_angles( Real & alpha, Real & beta, Real & gamma, Matrix M1, Matrix M2, bool const verbose = true );

void
create_euler_rotation(
	Matrix & M,
	Real const & alpha,
	Real const & beta,
	Real const & gamma,
	Vector const & /* axis1 not actually used*/,
	Vector const & axis2,
	Vector const & axis3
);

void
create_euler_rotation(
	Matrix & M,
	Real const & alpha,
	Real const & beta,
	Real const & gamma );

void
translate( pose::Pose & pose, Vector const & shift,
	pose::Pose const & ref_pose,
	utility::vector1< Size > const & moving_res );

void
rotate( pose::Pose & pose, Matrix const & M,
	pose::Pose const & ref_pose,
	utility::vector1< Size > const & moving_res,
	Vector const & centroid );

void
rotate( pose::Pose & pose, Matrix const & M,
	pose::Pose const & ref_pose,
	utility::vector1< Size > const & moving_res );

void
get_base_centroid_and_rotation_matrix( pose::Pose const & pose, Size const i, Vector & centroid, Matrix & M );

void
translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i, utility::vector1< Size > const & moving_res, bool const do_not_rotate = false );

void
translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i );


} //rigid_body
} //toolbox
} //protocols

#endif
