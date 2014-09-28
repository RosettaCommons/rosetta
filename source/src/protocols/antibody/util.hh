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
#include <core/pack/task/residue_selector/AndResidueSelector.hh>
#include <core/pack/task/residue_selector/NotResidueSelector.hh>
#include <core/pack/task/residue_selector/ChainSelector.hh>

//Protocol Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/clusters/CDRClusterEnumManager.hh>

//Utility Headers
#include <utility/vector1.hh>


using namespace core;
using namespace protocols::antibody::clusters;
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

///@brief Setup LH_A foldtree via docking.  Return dock_chains string.
std::string
setup_LH_A_foldtree(AntibodyInfoCOP ab_info, core::pose::Pose & pose);

///@brief Setup A_LH foldtree via docking. Return dock_chains string.
std::string
setup_A_LH_foldtree(AntibodyInfoCOP ab_info, core::pose::Pose & pose);


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

///@brief Get a set of loops for a boolean vector of CDRNameEnums including any stem residues.
protocols::loops::LoopsOP
get_cdr_loops(
	AntibodyInfoCOP ab_info,
	core::pose::Pose const & pose,
	utility::vector1<bool> cdrs,
	core::Size stem_size = 0 );


////////////////////////////////////////// things to compare to native //////////////////////////////////////////
core::Real
global_loop_rmsd ( const core::pose::Pose & pose_in,
                   const core::pose::Pose & native_pose,
                   loops::LoopsOP current_loop );


/// @brief return false if any cdr cutpoint is broken
bool
cutpoints_separation( core::pose::Pose & pose, AntibodyInfoOP & antibody_info );

/// @details Compute the separation at the cutpoint. The N-C distance of the
///   peptide bond which should be formed at the cutpoint. A closed loop is
///   assumed to have a gap < 1.9 Ang
core::Real
cutpoint_separation(core::pose::Pose & pose_in, Size cutpoint);

/////////////////////////////////// Epitope + Paratope ////////////////////////////////////////////////////////////////////

///@brief Get the epitope residues using the InterGroupNeighborsCalculator.  
vector1<bool>
select_epitope_residues(AntibodyInfoCOP ab_info, core::pose::Pose const & pose, core::Size const interface_distance = 10.0);

} //namespace antibody
} //namespace protocols


#endif //INCLUDED_protocols_antibody_util_HH




