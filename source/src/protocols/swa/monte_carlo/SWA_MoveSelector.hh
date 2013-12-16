// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/SWA_MoveSelector.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_monte_carlo_SWA_MoveSelector_HH
#define INCLUDED_protocols_swa_monte_carlo_SWA_MoveSelector_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/monte_carlo/SWA_MoveSelector.fwd.hh>
#include <protocols/swa/monte_carlo/SWA_Move.hh>
#include <core/pose/Pose.fwd.hh>

using namespace core;
using namespace core::pose;

namespace protocols {
namespace swa {
namespace monte_carlo {

	class SWA_MoveSelector: public utility::pointer::ReferenceCount {

	public:

		//constructor
		SWA_MoveSelector();

		//destructor
		~SWA_MoveSelector();

	public:

		void
		get_random_move_element_at_chain_terminus( pose::Pose & pose,
																							 SWA_Move & swa_move );

		void
		get_random_move_element_at_chain_terminus( pose::Pose & pose,
																							 SWA_Move & swa_move,
																							 utility::vector1< Size > const & sample_res /*leave empty if no filter*/);

		void
		get_resample_terminal_move_elements( pose::Pose & pose,
																				 utility::vector1< SWA_Move > & swa_moves );

		void
		get_resample_internal_move_elements( pose::Pose & pose,
																				 utility::vector1< SWA_Move > & swa_moves );

		void
		get_delete_move_elements( pose::Pose & pose,
															utility::vector1< SWA_Move > & swa_moves );

		Attachments
		get_attachments( pose::Pose & pose, Size const & moving_res );

		Attachments
		get_attachments( pose::Pose & pose, MoveElement const & move_element );

		void set_disallow_delete( bool const & setting ){ disallow_delete_ = setting; }
		bool disallow_delete() const{ return disallow_delete_; }

		void set_disallow_resample( bool const & setting ){ disallow_resample_ = setting; }
		bool disallow_resample() const{ return disallow_resample_; }

		void set_disallow_skip_bulge( bool const & setting ){ disallow_skip_bulge_ = setting; }
		bool disallow_skip_bulge() const{ return disallow_skip_bulge_; }

		void set_delete_terminal_only( bool const & setting ){ delete_terminal_only_ = setting; }
		bool delete_terminal_only() const{ return delete_terminal_only_; }

	private:

		utility::vector1< Size >
		get_partition_res( utility::vector1< bool > const & partition_definition, bool const setting );

		bool
		check_for_fixed_domain( pose::Pose const & pose,
														utility::vector1< Size> const & partition_res );

		void
		get_split_move_elements( pose::Pose & pose,
														 utility::vector1< SWA_Move > & swa_moves );

		void
		remove_from_consideration_first_multi_residue_move_element( utility::vector1< SWA_Move > & swa_moves,
																																bool remove_even_if_not_singlet );

		void
		get_terminal_move_elements( pose::Pose & pose,
																utility::vector1< SWA_Move > & swa_moves,
																MoveType const & move_type );

		void
		get_internal_move_elements( pose::Pose & pose,
																utility::vector1< SWA_Move > & swa_moves,
																MoveType const & move_type );

		void
		get_add_move_elements( pose::Pose & pose,
													 utility::vector1< SWA_Move > & swa_moves,
													 bool const disallow_skip_bulge );

		void
		filter_by_sample_res( utility::vector1< SWA_Move > & swa_moves,
													utility::vector1< Size > const & sample_res,
													utility::vector1< Size > const & res_list );

	private:

		bool disallow_delete_;
		bool disallow_resample_;
		bool disallow_skip_bulge_;
		bool delete_terminal_only_;

	};

} //monte_carlo
} //swa
} //protocols

#endif
