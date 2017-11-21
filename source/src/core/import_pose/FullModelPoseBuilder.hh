// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/FullModelPoseBuilder
/// @brief  create a pose from PDB in a way that sets up FullModelInfo
/// @author Andrew Watkins

#ifndef INCLUDED_core_import_pose_FullModelPoseBuilder_HH
#define INCLUDED_core_import_pose_FullModelPoseBuilder_HH

// Package headers
#include <core/import_pose/import_pose_options.fwd.hh>

// C++ headers
#include <iosfwd>

// Utility headers
#include <basic/Tracer.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileReaderOptions.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <tuple>

#include <utility/vector1.hh>

namespace core {
namespace import_pose {

class FullModelPoseBuilder {
public:
	FullModelPoseBuilder();
	~FullModelPoseBuilder() {}

	void set_options( utility::options::OptionCollection const & options ) { options_ = options; }
	void set_input_poses( utility::vector1< pose::PoseOP > input_poses ) { input_poses_ = input_poses; }
	void set_input_resnum_and_chain_and_segid( std::tuple< utility::vector1< Size >, utility::vector1< char >, utility::vector1< std::string > > const & input_resnum_and_chain_and_segid ) { input_resnum_and_chain_and_segid_ = input_resnum_and_chain_and_segid; }
	void set_cutpoint_open_in_full_model( utility::vector1< Size > const & cutpoint_open_in_full_model ) { cutpoint_open_in_full_model_ = cutpoint_open_in_full_model; }
	void set_fasta_file( std::string const & fasta_file ) { fasta_file_ = fasta_file; }
	void set_full_model_parameters( core::pose::full_model_info::FullModelParametersOP const & full_model_parameters ) { full_model_parameters_ = full_model_parameters; }

	void set_extra_minimize_res( utility::vector1< Size > const & extra_minimize_res ) { extra_minimize_res_ = extra_minimize_res; }
	void set_sample_res( utility::vector1< Size > const & sample_res ) { sample_res_ = sample_res; }
	void set_working_res( utility::vector1< Size > const & working_res ) { working_res_ = working_res; }
	void set_terminal_res( utility::vector1< Size > const & terminal_res ) { terminal_res_ = terminal_res; }
	void set_block_stack_above_res( utility::vector1< Size > const & block_stack_above_res ) { block_stack_above_res_ = block_stack_above_res; }
	void set_block_stack_below_res( utility::vector1< Size > const & block_stack_below_res ) { block_stack_below_res_ = block_stack_below_res; }
	void set_preferred_root_res( utility::vector1< Size > const & preferred_root_res ) { preferred_root_res_ = preferred_root_res; }
	void set_jump_res( utility::vector1< Size > const & jump_res ) { jump_res_ = jump_res; }
	void set_cutpoint_closed( utility::vector1< Size > const & cutpoint_closed ) { cutpoint_closed_ = cutpoint_closed; }
	void set_fiveprime_res( utility::vector1< Size > const & fiveprime_res ) { fiveprime_res_ = fiveprime_res; }
	void set_bulge_res( utility::vector1< Size > const & bulge_res ) { bulge_res_ = bulge_res; }
	void set_extra_minimize_jump_res( utility::vector1< Size > const & extra_minimize_jump_res );// { extra_minimize_jump_res_ = extra_minimize_jump_res; }
	void set_virtual_sugar_res( utility::vector1< Size > const & virtual_sugar_res ) { virtual_sugar_res_ = virtual_sugar_res; }
	void set_alignment_anchor_res( utility::vector1< Size > const & alignment_anchor_res ) { alignment_anchor_res_ = alignment_anchor_res; }
	void set_calc_rms_res( utility::vector1< Size > const & calc_rms_res ) { calc_rms_res_ = calc_rms_res; }
	void set_rna_syn_chi( utility::vector1< Size > const & rna_syn_chi ) { rna_syn_chi_ = rna_syn_chi; }
	void set_rna_anti_chi( utility::vector1< Size > const & rna_anti_chi ) { rna_anti_chi_ = rna_anti_chi; }
	void set_rna_north_sugar( utility::vector1< Size > const & rna_north_sugar ) { rna_north_sugar_ = rna_north_sugar; }
	void set_rna_south_sugar( utility::vector1< Size > const & rna_south_sugar ) { rna_south_sugar_ = rna_south_sugar; }
	void set_rna_sample_sugar( utility::vector1< Size > const & rna_sample_sugar ) { rna_sample_sugar_ = rna_sample_sugar; }
	void set_global_seq_file( std::string const & global_seq_file ) { global_seq_file_ = global_seq_file; }
	void set_disulfide_file( std::string const & disulfide_file ) { disulfide_file_ = disulfide_file; }
	void set_constraint_file( std::string const & constraint_file ) { constraint_file_ = constraint_file; }


	void initialize_input_poses_from_options( core::chemical::ResidueTypeSetCAP rsd_set );
	void initialize_further_from_options();
	void initialize_full_model_parameters();

	core::pose::PoseOP build();

private:

	void
	fill_full_model_info( core::pose::Pose & pose );
	void
	fill_full_model_info( utility::vector1< core::pose::PoseOP > & pose_ops );
	void
	fill_full_model_info( core::pose::Pose & pose, utility::vector1< core::pose::PoseOP > & other_pose_ops );
	void
	fill_full_model_info( utility::vector1< core::pose::Pose * > & pose_pointers );


	utility::options::OptionCollection options_;// = basic::options::option;
	utility::vector1< pose::PoseOP > input_poses_;
	std::tuple< utility::vector1< Size >, utility::vector1< char >, utility::vector1< std::string > > input_resnum_and_chain_and_segid_;
	utility::vector1< Size > cutpoint_open_in_full_model_;
	std::string fasta_file_ = "";
	core::pose::full_model_info::FullModelParametersOP full_model_parameters_ = nullptr;
	utility::vector1< Size > extra_minimize_res_;
	utility::vector1< Size > sample_res_;
	utility::vector1< Size > working_res_;
	utility::vector1< Size > terminal_res_;
	utility::vector1< Size > block_stack_above_res_;
	utility::vector1< Size > block_stack_below_res_;
	utility::vector1< Size > preferred_root_res_;
	utility::vector1< Size > jump_res_;
	utility::vector1< Size > cutpoint_closed_;
	utility::vector1< Size > fiveprime_res_;
	utility::vector1< Size > bulge_res_;
	utility::vector1< Size > extra_minimize_jump_res_;
	utility::vector1< Size > virtual_sugar_res_;
	utility::vector1< Size > alignment_anchor_res_;
	utility::vector1< Size > calc_rms_res_;
	utility::vector1< Size > rna_syn_chi_;
	utility::vector1< Size > rna_anti_chi_;
	utility::vector1< Size > rna_north_sugar_;
	utility::vector1< Size > rna_south_sugar_;
	utility::vector1< Size > rna_sample_sugar_;
	std::string global_seq_file_ = "";
	std::string disulfide_file_ = "";
	std::string constraint_file_ = "";

};

} // namespace import_pose
} // namespace core

#endif // INCLUDED_core_import_pose_FullModelPoseBuilder_HH
