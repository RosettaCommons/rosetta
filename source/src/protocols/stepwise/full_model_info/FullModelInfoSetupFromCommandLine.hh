// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/full_model_info/FullModelInfoSetupFromCommandLine.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_full_model_info_FullModelInfoSetupFromCommandLine_HH
#define INCLUDED_core_pose_full_model_info_FullModelInfoSetupFromCommandLine_HH

#include <protocols/stepwise/full_model_info/FullModelInfoSetupFromCommandLine.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>

using namespace core;
using utility::vector1;
using core::kinematics::FoldTree;

namespace core {
namespace pose {
namespace full_model_info {

	pose::PoseOP
	get_pdb_and_cleanup( std::string const input_file,
											 chemical::ResidueTypeSetCAP rsd_set );

	void
	get_other_poses( 	utility::vector1< pose::PoseOP > & other_poses,
										utility::vector1< std::string > const & other_files,
										chemical::ResidueTypeSetCAP rsd_set );

	pose::PoseOP
	initialize_pose_and_other_poses_from_command_line( core::chemical::ResidueTypeSetCAP rsd_set );

	void
	cleanup( pose::Pose & pose );

	void
	fill_full_model_info_from_command_line( pose::Pose & pose );

	void
	fill_full_model_info_from_command_line( utility::vector1< pose::PoseOP > & pose_ops );

	void
	fill_full_model_info_from_command_line( pose::Pose & pose, utility::vector1< PoseOP > & other_pose_ops );

	void
	fill_full_model_info_from_command_line( utility::vector1< pose::Pose * > & pose_pointers );

	void
	update_pose_fold_tree( pose::Pose & pose,
												 vector1< Size > const & res_list,
												 vector1< Size > const & extra_min_res,
												 vector1< Size > const & sample_res,
												 vector1< Size > const & jump_res,
												 vector1< Size > const & root_res );

	void
	define_chains( pose::Pose const & pose,
								 vector1< vector1< Size > > & all_res_in_chain,
								 vector1< vector1< Size > > & all_fixed_res_in_chain,
								 vector1< Size > const & res_list,
								 vector1< Size > const & extra_min_res );

	void
	setup_user_defined_jumps( vector1< Size > const & jump_res,
														vector1< Size > & jump_partners1,
														vector1< Size > & jump_partners2,
														vector1< vector1< bool > > & chain_connected,
														vector1< Size > const & res_list,
														vector1< vector1< Size > > const & all_res_in_chain );

	Size
	get_chain( Size const i, vector1< vector1< Size > > const & all_res_in_chain );

	void
	setup_jumps( vector1< Size > & jump_partners1,
							 vector1< Size > & jump_partners2,
							 vector1< vector1< bool > > & chain_connected,
							 vector1< vector1< Size > > const & all_res_in_chain );

	FoldTree
	get_tree( pose::Pose const & pose,
						vector1< Size > const & cuts,
						vector1< Size > const & jump_partners1,
						vector1< Size > const & jump_partners2 );

	FoldTree
	get_tree( Size const nres,
						vector1< Size > const & cuts,
						vector1< Size > const & jump_partners1,
						vector1< Size > const & jump_partners2,
						vector1< std::string > const & jump_atoms1,
						vector1< std::string > const & jump_atoms2 );

	void
	reroot( FoldTree & f, vector1< Size > const & res_list, vector1< Size > const & root_res );

	void
	add_cutpoint_closed( pose::Pose & pose,
											 vector1< Size > const & res_list,
											 vector1< Size > const & cutpoint_closed );

	void
	put_in_cutpoint( pose::Pose & pose, Size const i );

	void
	add_virtual_sugar_res( pose::Pose & pose,
												 vector1< Size > const & res_list,
												 vector1< Size > const & virtual_sugar_res );

} //full_model_info
} //pose
} //core

#endif
