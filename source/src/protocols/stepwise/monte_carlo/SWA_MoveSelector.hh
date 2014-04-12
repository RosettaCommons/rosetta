// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/SWA_MoveSelector.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_SWA_MoveSelector_HH
#define INCLUDED_protocols_stepwise_monte_carlo_SWA_MoveSelector_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/monte_carlo/SWA_MoveSelector.fwd.hh>
#include <protocols/stepwise/monte_carlo/SWA_Move.hh>
#include <core/pose/Pose.fwd.hh>

using namespace core;
using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace monte_carlo {

	class SWA_MoveSelector: public utility::pointer::ReferenceCount {

	public:

		//constructor
		SWA_MoveSelector();

		//destructor
		~SWA_MoveSelector();

	public:

		void
		get_random_add_or_delete_element( pose::Pose const & pose,
																			SWA_Move & swa_move );

		void
		get_random_add_or_delete_element( pose::Pose const & pose,
																			SWA_Move & swa_move,
																			utility::vector1< Size > const & sample_res /*leave empty if no filter*/);

		void
		get_resample_move_elements( pose::Pose const & pose,
																utility::vector1< SWA_Move > & swa_moves );

		void
		get_resample_terminal_move_elements( pose::Pose const & pose,
																				 utility::vector1< SWA_Move > & swa_moves );

		void
		get_resample_internal_move_elements( pose::Pose const & pose,
																				 utility::vector1< SWA_Move > & swa_moves );

		void
		get_resample_internal_local_move_elements( pose::Pose const & pose,
																							 utility::vector1< SWA_Move > & swa_moves );

		void
		get_intramolecular_delete_move_elements( pose::Pose const & pose,
																						 utility::vector1< SWA_Move > & swa_moves );

		Attachments
		get_attachments( pose::Pose const & pose, Size const & moving_res );

		Attachments
		get_attachments( pose::Pose const & pose, MoveElement const & move_element );

		void set_allow_delete( bool const & setting ){ allow_delete_ = setting; }
		bool allow_delete() const{ return allow_delete_; }

		void set_allow_skip_bulge( bool const & setting ){ allow_skip_bulge_ = setting; }
		bool allow_skip_bulge() const{ return allow_skip_bulge_; }

		void set_allow_internal_hinge( bool const & setting ){ allow_internal_hinge_ = setting; }
		bool allow_internal_hinge() const{ return allow_internal_hinge_; }

		void set_allow_internal_local( bool const & setting ){ allow_internal_local_ = setting; }
		bool allow_internal_local() const{ return allow_internal_local_; }

		void set_from_scratch_frequency( Real const & setting ){ from_scratch_frequency_ = setting; }
		Real from_scratch_frequency() const{ return from_scratch_frequency_; }

		void set_intermolecular_frequency( Real const & setting ){ intermolecular_frequency_ = setting; }
		Real intermolecular_frequency() const{ return intermolecular_frequency_; }

		void set_allow_shared_chains_in_dock_poses( bool const & setting ){ allow_shared_chains_in_dock_poses_ = setting; }
		bool allow_shared_chains_in_dock_poses() const{ return allow_shared_chains_in_dock_poses_; }

	private:

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void
		get_intramolecular_add_and_delete_elements( pose::Pose const & pose,
																								utility::vector1< SWA_Move > & swa_moves,
																								utility::vector1< Size > const & sample_res /*leave empty if no filter*/);

		void
		get_intramolecular_add_move_elements( pose::Pose const & pose,
																					utility::vector1< SWA_Move > & swa_moves );

		void
		get_intramolecular_split_move_elements( pose::Pose const & pose,
																						utility::vector1< SWA_Move > & swa_moves,
																						MoveType const move_type );

		void
		remove_from_consideration_first_multi_residue_move_element( utility::vector1< SWA_Move > & swa_moves,
																																bool remove_even_if_not_singlet );

		void
		filter_by_sample_res( utility::vector1< SWA_Move > & swa_moves,
													utility::vector1< Size > const & sample_res );

		bool
		check_for_fixed_domain_or_from_scratch(  pose::Pose const & pose,
																						 utility::vector1< Size> const & partition_res ) const;

		bool
		already_instantiated_in_pose( pose::Pose const & pose, Size const & resnum_in_full_model_numbering ) const;

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void
		get_intermolecular_add_and_delete_elements( pose::Pose const & pose,
																								utility::vector1< SWA_Move > & swa_moves,
																								utility::vector1< Size > const & sample_res /*leave empty if no filter*/);

		void
		get_intermolecular_add_move_elements( pose::Pose const & pose,
																					utility::vector1< SWA_Move > & swa_moves );

		void
		get_intermolecular_delete_move_elements( pose::Pose const & pose,
																						utility::vector1< SWA_Move > & swa_moves );

		void
		get_intermolecular_split_move_elements( pose::Pose const & pose,
																						utility::vector1< SWA_Move > & swa_moves,
																						MoveType const move_type );

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void
		get_from_scratch_add_and_delete_elements( pose::Pose const & pose,
																									 utility::vector1< SWA_Move > & swa_moves,
																									 utility::vector1< Size > const & sample_res /*leave empty if no filter*/);

		void
		get_from_scratch_add_move_elements( pose::Pose const & pose,
																				utility::vector1< SWA_Move > & swa_moves );

		void
		figure_out_from_scratch_delete_move_elements( pose::Pose const & pose,
																									utility::vector1< SWA_Move > & swa_moves ) const;

		void
		get_from_scratch_delete_move_elements( pose::Pose const & pose,
																					 utility::vector1< SWA_Move > & swa_moves );

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// other moves [used in resampling]
		void
		get_terminal_move_elements( pose::Pose const & pose,
																utility::vector1< SWA_Move > & swa_moves,
																MoveType const & move_type );

		void
		get_internal_move_elements( pose::Pose const & pose,
																utility::vector1< SWA_Move > & swa_moves,
																MoveType const & move_type );


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// helper functions
		bool
		partition_splits_a_fixed_domain( utility::vector1< Size > const & partition_definition,
																		 utility::vector1< Size > const & domain_map ) const;

		utility::vector1< Size >
		get_partition_res( utility::vector1< bool > const & partition_definition, bool const setting ) const;

		utility::vector1< Size >
		get_unique_chains( pose::Pose const & pose );

		bool
		share_chains( utility::vector1< Size > const & current_unique_chains,
									pose::Pose const & pose,
									Size const j_full );


	private:

		bool allow_delete_;
		bool allow_skip_bulge_;
		bool allow_internal_hinge_;
		bool allow_internal_local_;
		Real from_scratch_frequency_;
		Real intermolecular_frequency_;
		bool only_dock_preexisting_chunks_;
		bool allow_shared_chains_in_dock_poses_;
		bool resampling_;
	};

} //monte_carlo
} //stepwise
} //protocols

#endif
