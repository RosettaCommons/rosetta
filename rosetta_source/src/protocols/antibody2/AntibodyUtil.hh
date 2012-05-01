// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody2/AntibodyUtil.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#ifndef INCLUDED_protocols_antibody2_AntibodyUtil_hh
#define INCLUDED_protocols_antibody2_AntibodyUtil_hh


#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/antibody2/AntibodyInfo.hh>



using namespace core;
///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace antibody2 {
void simple_one_loop_fold_tree(
                                   core::pose::Pose & pose,
                                   loops::Loop const & loop);
    
    
    
    
void simple_fold_tree(
                          core::pose::Pose & pose_in,
                          core::Size jumppoint1,
                          core::Size cutpoint,
                          core::Size jumppoint2);
    
    
    
void setup_simple_fold_tree(
                                core::Size jumppoint1,
                                core::Size cutpoint,
                                core::Size jumppoint2,
                                core::Size nres,
                                core::pose::Pose & pose_in );
    
    
void all_cdr_VL_VH_fold_tree( core::pose::Pose & pose_in, 
                            const loops::Loops & loops );
    
    
bool CDR_H3_filter(
                       const core::pose::Pose & pose_in,
                       loops::Loop & input_loop,
                       bool is_camelid);
    
    
core::pack::task::TaskFactoryOP setup_packer_task( core::pose::Pose & pose_in);

/*    void dle_extreme_repack(pose::Pose & pose_in,
                            int repack_cycles,
                            ObjexxFCL::FArray1D_bool & allow_repack,
                            bool rt_min,
                            bool rotamer_trials,
                            bool force_one_repack,
                            bool use_unbounds);
  */  

/// @brief return false if any cdr cutpoint is broken
bool cutpoints_separation( core::pose::Pose & pose, AntibodyInfoOP & antibody_info );
    
    
    
    
// Compute the separation at the cutpoint. The N-C distance of the
// peptide bond which should be formed at the cutpoint. A closed loop is
// assumed to have a gap < 1.9 Ang
core::Real cutpoint_separation(core::pose::Pose & pose_in, Size cutpoint);

    
    
    
core::Real global_loop_rmsd ( const core::pose::Pose & pose_in, 
                              const core::pose::Pose & native_pose, 
                              loops::LoopOP current_loop );

    
std::string get_seq_from_a_loop(core::pose::Pose & pose_in, loops::LoopOP  loop);
    

} //namespace antibody2
} //namespace protocols


#endif //INCLUDED_protocols_loops_AntibodyUtil_HH




