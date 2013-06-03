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
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_AntibodyUtil_hh
#define INCLUDED_protocols_antibody_AntibodyUtil_hh


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
    
/// @brief calculates the VH_VL packing angle from 2 sheet definitions on each chain from ab_info.
core::Real
vl_vh_packing_angle ( const core::pose::Pose & pose_in, AntibodyInfoOP ab_info );

/// @brief return false if any cdr cutpoint is broken
bool
cutpoints_separation( core::pose::Pose & pose, AntibodyInfoOP & antibody_info );
    
/// @details Compute the separation at the cutpoint. The N-C distance of the
///   peptide bond which should be formed at the cutpoint. A closed loop is
///   assumed to have a gap < 1.9 Ang
core::Real
cutpoint_separation(core::pose::Pose & pose_in, Size cutpoint);


///// kink measures /////

// @brief returns distance of the sc-sc Hbond across the strands at the beginning of the H3 kink (typically Asp-Arg)
core::Real
kink_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & abinfo);

// @brief returns distance for the bb-bb Hbond across the strands at the begining of the kink (typically Asp-Arg)
core::Real
kink_bb_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & abinfo);


// @brief returns distance of the Trp sc-bb Hbond across the H3 kink residues (n-1 to n+2)
core::Real
kink_Trp_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & abinfo);

// @brief returns a pair of reals for the q distance and qbond dihedral angle from the four kink residues of the H3 C-terminal end
std::pair<core::Real,core::Real>
kink_dihedral( const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & abinfo, bool debug=false);



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CDR Clusters
//
//

/// @brief Sets dihedral harmonic constraints to CDR and scorefxn using info in AntibodyInfo
/// @details Currently requires Modified_AHO numbering.  A server will be available soon. 
void
set_harmonic_constraints(AntibodyInfoOP & ab_info, pose::Pose & pose, core::scoring::ScoreFunctionOP & scorefxn);

/// @brief Sets dihedral harmonic constraints to CDR using cluster info in AntibodyInfo
/// @details Currently requires Modified_AHO numbering.  A server will be available soon. 
void
set_harmonic_constraints(AntibodyInfoOP & ab_info, pose::Pose & pose);

/// @brief set a harmonic constraint to a CDR based on cluster type
/// @details Currently requires Modified_AHO numbering.  A server will be available soon. 
void
set_harmonic_constraint(AntibodyInfoOP & ab_info, pose::Pose & pose, CDRClusterEnum const cluster);
	

	

	
} //namespace antibody
} //namespace protocols


#endif //INCLUDED_protocols_antibody_AntibodyUtil_HH




