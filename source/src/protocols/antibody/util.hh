// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/util.hh
/// @brief Utility functions for the Antibody namespace
/// @author Jianqing Xu (xubest@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_util_hh
#define INCLUDED_protocols_antibody_util_hh


//Core Headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

//Protocol Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/CDRClusterEnumManager.hh>

//Utility Headers
#include <utility/vector1.hh>


using namespace core;
///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace antibody {

void
simple_one_loop_fold_tree(
    core::pose::Pose & pose,
    loops::Loop const & loop);

void
simple_fold_tree(
    core::pose::Pose & pose_in,
    core::Size jumppoint1,
    core::Size cutpoint,
    core::Size jumppoint2);


/// @brief Very specific packertask,
/// @details Ubound rotamer options, repack only, protein only, no disulfides.
core::pack::task::TaskFactoryOP
setup_packer_task( core::pose::Pose & pose_in);


/// @brief align current Fv to native.Fv
void
align_to_native( core::pose::Pose & pose,
                 core::pose::Pose const & native_pose,
                 AntibodyInfoOP const ab_info,
                 AntibodyInfoOP const native_ab_info,
                 std::string const & request_chain="LH");

bool
CDR_H3_filter_legacy_code_with_old_rule(
    const core::pose::Pose & pose_in,
    loops::Loop & input_loop,
    bool is_camelid);

bool
CDR_H3_cter_filter(
    const core::pose::Pose & pose_in,
    AntibodyInfoOP ab_info);


////////////////////////////////////////// things to compare to native //////////////////////////////////////////
core::Real
global_loop_rmsd ( const core::pose::Pose & pose_in,
                   const core::pose::Pose & native_pose,
                   loops::LoopsOP current_loop );


////////////////////////////////////////// antibody metrics //////////////////////////////////////////

/// @brief return false if any cdr cutpoint is broken
bool
cutpoints_separation( core::pose::Pose & pose, AntibodyInfoOP & antibody_info );

/// @details Compute the separation at the cutpoint. The N-C distance of the
///   peptide bond which should be formed at the cutpoint. A closed loop is
///   assumed to have a gap < 1.9 Ang
core::Real
cutpoint_separation(core::pose::Pose & pose_in, Size cutpoint);



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CDR Clusters
//
//


/// @brief Adds dihedral harmonic constraints to Pose CDRs using cluster info in AntibodyInfo
/// @details Currently requires Modified_AHO numbering. Returns map of success/failure
std::map<CDRNameEnum, bool>
add_harmonic_cluster_constraints(AntibodyInfoOP ab_info, pose::Pose & pose);

///@brief Same as above, but adds constraints to the vector so they can be identified and removed from the pose if needed.
///@details Returns map of success/failure
std::map<CDRNameEnum, bool>
add_harmonic_cluster_constraints(AntibodyInfoOP ab_info, pose::Pose & pose, utility::vector1< core::scoring::constraints::ConstraintCOP > constraints);


/// @brief Adds a harmonic constraint to a Pose CDR based on cluster type
/// @details Currently requires Modified_AHO numbering.
bool
add_harmonic_cluster_constraint(AntibodyInfoCOP ab_info, pose::Pose & pose, CDRClusterEnum const cluster);

///@brief Same as above, but adds constraints to the vector so they can be identified and removed from the pose if needed.
///@details Returns true or false depending on success
bool
add_harmonic_cluster_constraint(AntibodyInfoCOP ab_info, pose::Pose & pose, CDRClusterEnum const cluster, utility::vector1< core::scoring::constraints::ConstraintCOP > constraints);

/// @brief Uses ConstraintSetMover to set a constraint.  Through cstmover, only 1 constraint can be set a t a time.
/// @details Currently requires Modified_AHO numbering.  Returns True or False depending on success.
bool
set_harmonic_cluster_constraint(AntibodyInfoCOP ab_info, pose::Pose & pose, CDRClusterEnum const cluster);


/// @brief Gets the cluster constraint name.  Returns NA if not found.
std::string
get_harmonic_cluster_constraint_filename(AntibodyInfoCOP ab_info, CDRClusterEnum const cluster);


///@brief Very basic way to check to make sure pose residues are Modified_AHO (North, et al) scheme, which allows the clustering.
///@details If any of these anchor residues that are checked are missing, it will return false.
bool
check_if_pose_renumbered_for_clusters(pose::Pose const & pose);
	


} //namespace antibody
} //namespace protocols


#endif //INCLUDED_protocols_antibody_util_HH




