// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/AntibodyUtil.hh
/// @brief Utility functions for the Antibody namespace
/// @author Jianqing Xu (xubest@gmail.com)

#ifndef INCLUDED_protocols_antibody_AntibodyUtil_hh
#define INCLUDED_protocols_antibody_AntibodyUtil_hh


#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/antibody/AntibodyInfo.hh>



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

/// @brief align current Fv to native.Fv
void
align_to_native( core::pose::Pose & pose, 
	core::pose::Pose & native_pose,
	AntibodyInfoOP ab_info,
	AntibodyInfoOP native_ab_info);

bool 
CDR_H3_filter_legacy_code_with_old_rule(
	const core::pose::Pose & pose_in,
	loops::Loop & input_loop,
	bool is_camelid);
    
bool 
CDR_H3_cter_filter(
	const core::pose::Pose & pose_in,
	AntibodyInfoOP ab_info);

core::Real
global_loop_rmsd ( const core::pose::Pose & pose_in, 
	const core::pose::Pose & native_pose, 
	loops::LoopsOP current_loop );

/// @brief calculates the VH_VL packing angle from 2 sheet definitions on each chain from ab_info.
core::Real
vl_vh_packing_angle ( const core::pose::Pose & pose_in, AntibodyInfoOP ab_info );

/// @brief Very specific packertask,
/// @details Ubound rotamer options, repack only, protein only, no disulfides.
core::pack::task::TaskFactoryOP
setup_packer_task( core::pose::Pose & pose_in);

/// @brief return false if any cdr cutpoint is broken
bool
cutpoints_separation( core::pose::Pose & pose, AntibodyInfoOP & antibody_info );
    
/// @details Compute the separation at the cutpoint. The N-C distance of the
///   peptide bond which should be formed at the cutpoint. A closed loop is
///   assumed to have a gap < 1.9 Ang
core::Real
cutpoint_separation(core::pose::Pose & pose_in, Size cutpoint);

/*    void dle_extreme_repack(pose::Pose & pose_in,
	int repack_cycles,
	ObjexxFCL::FArray1D_bool & allow_repack,
	bool rt_min,
	bool rotamer_trials,
	bool force_one_repack,
	bool use_unbounds);
  */  

} //namespace antibody
} //namespace protocols


#endif //INCLUDED_protocols_antibody_AntibodyUtil_HH




