// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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


namespace protocols {
namespace toolbox {
namespace rigid_body {

typedef  numeric::xyzMatrix< core::Real > Matrix;

void
get_atom_coordinates( utility::vector1< std::pair < core::id::AtomID, numeric::xyzVector< core::Real > > > & xyz_list,
	core::Size const & seq_num,
	core::conformation::Residue const & rsd_at_origin,
	core::kinematics::Stub const & moving_res_base_stub );

void
get_specific_atom_coordinate( core::Vector & atom_pos,
	core::kinematics::Stub const & moving_res_base_stub );

void
get_specific_atom_coordinate( core::Size const atom_index,
	core::Vector & atom_pos,
	core::conformation::Residue const & rsd_at_origin,
	core::kinematics::Stub const & moving_res_base_stub );


void
get_specific_atom_coordinate( std::string const & atom_name,
	core::Vector & atom_pos,
	core::conformation::Residue const & rsd_at_origin,
	core::kinematics::Stub const & moving_res_base_stub );


/////////////////////////////////////////////////////////////
// following functions are helpers for RigidBodyStepWiseSampler
core::Size
figure_out_reference_res_for_jump( core::pose::Pose const & pose,
	core::Size const moving_res );

utility::vector1< core::Size >
figure_out_moving_partition_res( core::pose::Pose const & pose,
	core::Size const moving_res, core::Size const reference_res );

void
set_to_origin( core::pose::Pose & pose, core::Size const seq_num, bool verbose );

core::kinematics::Stub
initialize_stub( core::conformation::Residue const & rsd );

core::pose::PoseCOP
transform_moving_partition_to_origin( core::pose::Pose const & pose_start,
	core::Size const moving_res,
	utility::vector1< core::Size > const & moving_partition_res_ );

/////////////////////////////////////////////////////////////////////////////////////////
// Following were developed for a rigid-body rotation enumeration framework for
// computing partition functions for the RNA nearest-neighbor rules. No
// longer in use. However, could be very useful in lo-res sampling, so I didn't
// remove. If these are not in use by 2016, delete!
void
get_euler_angles( core::Real & alpha, core::Real & beta, core::Real & gamma, Matrix M1, Matrix M2, bool const verbose = true );

void
create_euler_rotation(
	Matrix & M,
	core::Real const & alpha,
	core::Real const & beta,
	core::Real const & gamma,
	core::Vector const & /* axis1 not actually used*/,
	core::Vector const & axis2,
	core::Vector const & axis3
);

void
create_euler_rotation(
	Matrix & M,
	core::Real const & alpha,
	core::Real const & beta,
	core::Real const & gamma );

void
translate( core::pose::Pose & pose, core::Vector const & shift,
	core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & moving_res );

void
rotate( core::pose::Pose & pose, Matrix const & M,
	core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & moving_res,
	core::Vector const & centroid );

void
rotate( core::pose::Pose & pose, Matrix const & M,
	core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & moving_res );

void
get_base_centroid_and_rotation_matrix( core::pose::Pose const & pose, core::Size const i, core::Vector & centroid, Matrix & M );

void
translate_and_rotate_residue_to_origin( core::pose::Pose & pose, core::Size const i, utility::vector1< core::Size > const & moving_res, bool const do_not_rotate = false );

void
translate_and_rotate_residue_to_origin( core::pose::Pose & pose, core::Size const i );


} //rigid_body
} //toolbox
} //protocols

#endif
