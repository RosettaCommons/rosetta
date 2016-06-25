// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace core::pose;

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
	figure_out_all_possible_moves( pose::Pose const & pose, bool const verbose = false );

	void
	output_moves() const;

	StepWiseMove
	select_random_move( pose::Pose const & pose ) const;

	utility::vector1< StepWiseMove > const & swa_moves() const { return swa_moves_; }

	Real
	proposal_probability( StepWiseMove const & swa_move ) const;

	void
	get_add_or_delete_element( pose::Pose const & pose,
		StepWiseMove & swa_move );

	void
	get_resample_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		bool const save_moves = true );

	void
	get_resample_terminal_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_resample_internal_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_resample_internal_local_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_intramolecular_delete_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	Attachments
	get_attachments( pose::Pose const & pose, Size const & moving_res );

	Attachments
	get_attachments( pose::Pose const & pose, MoveElement const & move_element );

	void set_allow_delete( bool const & setting ){ allow_delete_ = setting; }
	bool allow_delete() const{ return allow_delete_; }

	void set_choose_random( bool const & setting ){ choose_random_ = setting; }
	bool choose_random() const{ return choose_random_; }

	void set_force_unique_moves( bool const & setting ){ force_unique_moves_ = setting; }
	bool force_unique_moves() const{ return force_unique_moves_; }

	monte_carlo::submotif::SubMotifLibraryCOP submotif_library() { return submotif_library_; }
	void set_submotif_library( monte_carlo::submotif::SubMotifLibraryCOP setting ) { submotif_library_ = setting; }

	StepWiseMove
	reverse_move( StepWiseMove const & swa_move, pose::Pose const & pose_before, pose::Pose const & pose_after ) const;

	bool
	just_simple_cycles( StepWiseMove const & swa_move, pose::Pose const & pose,
		bool const verbose = false ) const;

	void
	set_options( protocols::stepwise::monte_carlo::mover::options::StepWiseMoveSelectorOptionsCOP setting ){ options_ = setting; }
	protocols::stepwise::monte_carlo::mover::options::StepWiseMoveSelectorOptionsCOP options() const { return options_; }


private:

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	fill_moves_for_pose( pose::Pose const & pose );

	void
	fill_moves_for_other_poses( pose::Pose const & pose );

	void
	fill_denovo_moves( pose::Pose const & pose );

	void
	fill_vary_loop_length_moves( pose::Pose const & pose );

	void
	save_moves( utility::vector1< StepWiseMove > const & moves,
		Real const total_weight = 1,
		bool const clear_moves_before_saving = false );

	Real
	sum_probabilities( Size const start_idx = 1,
		Size const final_idx = 0 ) const;

	void
	normalize_probabilities( Size const start_idx = 1,
		Size const final_idx = 0,
		Real const desired_weight = 1 );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_intramolecular_add_and_delete_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_intramolecular_add_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_intramolecular_split_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		MoveType const move_type );

	bool
	check_from_scratch(  pose::Pose const & pose,
		utility::vector1< bool > const & partition_definition) const;

	void
	remove_from_consideration_first_multi_residue_move_element( utility::vector1< StepWiseMove > & swa_moves,
		bool remove_even_if_not_singlet );

	bool
	is_addable_res( Size const n, pose::Pose const & pose ) const;

	bool
	remnant_would_be_deleted(
		pose::Pose const & pose,
		utility::vector1 < Size > const & partition ) const;

	bool
	both_remnants_would_be_deleted(
		pose::Pose const & pose,
		utility::vector1 < Size > const & partition1,
		utility::vector1 < Size > const & partition2 ) const;

	bool
	partitions_split_a_submotif( pose::Pose const & pose, utility::vector1< Size > const & partition1, utility::vector1< Size > const & partition2 ) const;

	bool
	check_for_intramolecular_submotif_jump( pose::Pose const & pose, Size const & moving_res, Size const & attached_res ) const;

	bool
	check_for_input_domain_or_from_scratch(  pose::Pose const & pose,
		utility::vector1< Size> const & partition_res ) const;

	bool
	already_instantiated_in_pose( pose::Pose const & pose, Size const & resnum_in_full_model_numbering ) const;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_docking_add_and_delete_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	std::set< Size >
	get_dock_domains( utility::vector1< Size > const & move_element,
		pose::Pose const & pose ) const;

	void
	figure_out_already_docked(
		pose::Pose const & pose,
		std::map< std::pair< Size, Size >, bool > & already_docked ) const;

	void
	get_docking_add_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_docking_delete_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_docking_split_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		MoveType const move_type );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_submotif_add_moves( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	filter_add_submotif_moves_to_not_redock_domain( utility::vector1< StepWiseMove > & submotif_add_moves, pose::Pose const & pose ) const;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_from_scratch_add_and_delete_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_from_scratch_add_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	figure_out_from_scratch_delete_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves ) const;

	void
	get_from_scratch_delete_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves );

	void
	get_vary_loop_length_moves( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & vary_loop_length_moves ) const;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// other moves [used in remodeler]
	void
	get_terminal_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		MoveType const & move_type );

	void
	get_internal_move_elements( pose::Pose const & pose,
		utility::vector1< StepWiseMove > & swa_moves,
		MoveType const & move_type );


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// helper functions
	bool
	partition_splits_an_input_domain( utility::vector1< Size > const & partition_definition,
		utility::vector1< Size > const & domain_map ) const;

	utility::vector1< Size >
	get_partition_res( utility::vector1< bool > const & partition_definition, bool const setting ) const;

	utility::vector1< Size >
	get_unique_chains( pose::Pose const & pose );

	bool
	share_chains( utility::vector1< Size > const & current_unique_chains,
		pose::Pose const & pose,
		Size const j_full );

	utility::vector1< Size >
	figure_out_first_res_in_pose( utility::vector1< Size > const & pose_domain_map ) const;

	void
	filter_pose_order( utility::vector1< StepWiseMove > & swa_moves,
		pose::Pose const & pose ) const;


	//////////////////////////////////////////////////////////////////////////////////////////
	// reverse_moves
	Size
	get_actual_moving_res( StepWiseMove const & swa_move, pose::Pose const & pose ) const;

	StepWiseMove
	reverse_add_move( StepWiseMove const & swa_move, pose::Pose const & pose_after ) const;

	StepWiseMove
	reverse_add_submotif_move( StepWiseMove const & swa_move, pose::Pose const & pose_after ) const;

	StepWiseMove
	reverse_resample_move( StepWiseMove const & swa_move, pose::Pose const & pose_after ) const;

	StepWiseMove
	reverse_delete_move( StepWiseMove const & swa_move, pose::Pose const & pose_before, pose::Pose const & pose_after ) const;

	Attachment
	figure_out_attachment( Size const moving_res, Size const attached_res, pose::Pose const & pose ) const;

	StepWiseMove
	ordered_move_from_partition( Size const actual_moving_res, Size const actual_attached_res,
		pose::Pose const & pose, MoveType const & move_type ) const;

	void
	filter_complex_cycles( utility::vector1< StepWiseMove > & swa_moves, pose::Pose const & pose ) const;

private:

	protocols::stepwise::monte_carlo::mover::options::StepWiseMoveSelectorOptionsCOP options_;

	bool allow_delete_;
	bool force_unique_moves_;
	bool choose_random_;

	monte_carlo::submotif::SubMotifLibraryCOP submotif_library_;

	utility::vector1< StepWiseMove > swa_moves_;
	utility::vector1< Real > proposal_probabilities_;

};

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
