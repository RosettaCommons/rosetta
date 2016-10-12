// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_StepWiseMoveSelector_HH
#define INCLUDED_protocols_stepwise_monte_carlo_StepWiseMoveSelector_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/options/StepWiseMoveSelectorOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <set>
#include <map>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

class StepWiseMoveSelector: public utility::pointer::ReferenceCount {

public:

	//constructor
	StepWiseMoveSelector();

	//constructor
	StepWiseMoveSelector( options::StepWiseMoveSelectorOptionsCOP options );

	//constructor
	StepWiseMoveSelector( StepWiseMoveSelector const & src );

	//destructor
	~StepWiseMoveSelector();

public:

	StepWiseMoveSelectorOP clone() const;

	bool
	figure_out_all_possible_moves( core::pose::Pose const & pose, bool const verbose = false );

	void
	output_moves() const;

	StepWiseMove
	select_random_move( core::pose::Pose const & pose ) const;

	utility::vector1< StepWiseMove > const & swa_moves() const { return swa_moves_; }

	core::Real
	proposal_probability( StepWiseMove const & swa_move ) const;

	void
	get_add_or_delete_element( core::pose::Pose const & pose,
		StepWiseMove & swa_move );

	void
	get_resample_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		bool const save_moves = true );

	void
	get_resample_terminal_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_resample_internal_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_resample_internal_local_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_intramolecular_delete_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	Attachments
	get_attachments( core::pose::Pose const & pose, core::Size const & moving_res );

	Attachments
	get_attachments( core::pose::Pose const & pose, MoveElement const & move_element );

	void set_allow_delete( bool const & setting ){ allow_delete_ = setting; }
	bool allow_delete() const{ return allow_delete_; }

	void set_choose_random( bool const & setting ){ choose_random_ = setting; }
	bool choose_random() const{ return choose_random_; }

	void set_force_unique_moves( bool const & setting ){ force_unique_moves_ = setting; }
	bool force_unique_moves() const{ return force_unique_moves_; }

	monte_carlo::submotif::SubMotifLibraryCOP submotif_library() { return submotif_library_; }
	void set_submotif_library( monte_carlo::submotif::SubMotifLibraryCOP setting ) { submotif_library_ = setting; }

	StepWiseMove
	reverse_move( StepWiseMove const & swa_move, core::pose::Pose const & pose_before, core::pose::Pose const & pose_after ) const;

	bool
	just_simple_cycles( StepWiseMove const & swa_move, core::pose::Pose const & pose,
		bool const verbose = false ) const;

	void
	set_options( protocols::stepwise::monte_carlo::mover::options::StepWiseMoveSelectorOptionsCOP setting ){ options_ = setting; }
	protocols::stepwise::monte_carlo::mover::options::StepWiseMoveSelectorOptionsCOP options() const { return options_; }


private:

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	fill_moves_for_pose( core::pose::Pose const & pose );

	void
	fill_moves_for_other_poses( core::pose::Pose const & pose );

	void
	fill_denovo_moves( core::pose::Pose const & pose );

	void
	fill_vary_loop_length_moves( core::pose::Pose const & pose );

	void
	save_moves( utility::vector1< StepWiseMove > const & moves,
		core::Real const total_weight = 1,
		bool const clear_moves_before_saving = false );

	core::Real
	sum_probabilities( core::Size const start_idx = 1,
		core::Size const final_idx = 0 ) const;

	void
	normalize_probabilities( core::Size const start_idx = 1,
		core::Size const final_idx = 0,
		core::Real const desired_weight = 1 );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_intramolecular_add_and_delete_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_intramolecular_add_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_intramolecular_split_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		MoveType const move_type );

	bool
	check_from_scratch(  core::pose::Pose const & pose,
		utility::vector1< bool > const & partition_definition) const;

	void
	remove_from_consideration_first_multi_residue_move_element( utility::vector1< StepWiseMove > & swa_moves,
		bool remove_even_if_not_singlet );

	bool
	is_addable_res( core::Size const n, core::pose::Pose const & pose ) const;

	bool
	remnant_would_be_deleted(
		core::pose::Pose const & pose,
		utility::vector1 < core::Size > const & partition ) const;

	bool
	both_remnants_would_be_deleted(
		core::pose::Pose const & pose,
		utility::vector1 < core::Size > const & partition1,
		utility::vector1 < core::Size > const & partition2 ) const;

	bool
	partitions_split_a_submotif( core::pose::Pose const & pose, utility::vector1< core::Size > const & partition1, utility::vector1< core::Size > const & partition2 ) const;

	bool
	check_for_intramolecular_submotif_jump( core::pose::Pose const & pose, core::Size const & moving_res, core::Size const & attached_res ) const;

	bool
	check_for_input_domain_or_from_scratch(  core::pose::Pose const & pose,
		utility::vector1< core::Size> const & partition_res ) const;

	bool
	already_instantiated_in_pose( core::pose::Pose const & pose, core::Size const & resnum_in_full_model_numbering ) const;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_docking_add_and_delete_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	std::set< core::Size >
	get_dock_domains( utility::vector1< core::Size > const & move_element,
		core::pose::Pose const & pose ) const;

	void
	figure_out_already_docked(
		core::pose::Pose const & pose,
		std::map< std::pair< core::Size, core::Size >, bool > & already_docked ) const;

	void
	get_docking_add_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_docking_delete_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_docking_split_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		MoveType const move_type );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_submotif_add_moves( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	filter_add_submotif_moves_to_not_redock_domain( utility::vector1< StepWiseMove > & submotif_add_moves, core::pose::Pose const & pose ) const;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_from_scratch_add_and_delete_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_from_scratch_add_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	figure_out_from_scratch_delete_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves ) const;

	void
	get_from_scratch_delete_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_vary_loop_length_moves( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & vary_loop_length_moves ) const;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// other moves [used in remodeler]
	void
	get_terminal_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		MoveType const & move_type );

	void
	get_internal_move_elements( core::pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		MoveType const & move_type );


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// helper functions
	bool
	partition_splits_an_input_domain( utility::vector1< core::Size > const & partition_definition,
		utility::vector1< core::Size > const & domain_map ) const;

	utility::vector1< core::Size >
	get_partition_res( utility::vector1< bool > const & partition_definition, bool const setting ) const;

	utility::vector1< core::Size >
	get_unique_chains( core::pose::Pose const & pose );

	bool
	share_chains( utility::vector1< core::Size > const & current_unique_chains,
		core::pose::Pose const & pose,
		core::Size const j_full );

	utility::vector1< core::Size >
	figure_out_first_res_in_pose( utility::vector1< core::Size > const & pose_domain_map ) const;

	void
	filter_pose_order( utility::vector1< StepWiseMove > & swa_moves,
		core::pose::Pose const & pose ) const;


	//////////////////////////////////////////////////////////////////////////////////////////
	// reverse_moves
	core::Size
	get_actual_moving_res( StepWiseMove const & swa_move, core::pose::Pose const & pose ) const;

	StepWiseMove
	reverse_add_move( StepWiseMove const & swa_move, core::pose::Pose const & pose_after ) const;

	StepWiseMove
	reverse_add_submotif_move( StepWiseMove const & swa_move, core::pose::Pose const & pose_after ) const;

	StepWiseMove
	reverse_resample_move( StepWiseMove const & swa_move, core::pose::Pose const & pose_after ) const;

	StepWiseMove
	reverse_delete_move( StepWiseMove const & swa_move, core::pose::Pose const & pose_before, core::pose::Pose const & pose_after ) const;

	Attachment
	figure_out_attachment( core::Size const moving_res, core::Size const attached_res, core::pose::Pose const & pose ) const;

	StepWiseMove
	ordered_move_from_partition( core::Size const actual_moving_res, core::Size const actual_attached_res,
		core::pose::Pose const & pose, MoveType const & move_type ) const;

	void
	filter_complex_cycles( utility::vector1< StepWiseMove > & swa_moves, core::pose::Pose const & pose ) const;

private:

	protocols::stepwise::monte_carlo::mover::options::StepWiseMoveSelectorOptionsCOP options_;

	bool allow_delete_;
	bool force_unique_moves_;
	bool choose_random_;

	monte_carlo::submotif::SubMotifLibraryCOP submotif_library_;

	utility::vector1< StepWiseMove > swa_moves_;
	utility::vector1< core::Real > proposal_probabilities_;

};

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
