// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_MonteCarloUtil.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_monte_carlo_SWA_MonteCarloUtil_hh
#define INCLUDED_protocols_swa_monte_carlo_SWA_MonteCarloUtil_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/swa/monte_carlo/SWA_Move.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace monte_carlo {

	MovingResidueCase
	get_moving_residue_case( pose::Pose const & pose, Size const i );

	MovingResidueCase
	get_moving_residue_case( pose::Pose const & pose, utility::vector1< Size > const & moving_chunk );

	void
	get_potential_delete_chunks( pose::Pose & pose,
															 utility::vector1< SWA_Move > & swa_moves );

	void
	get_potential_resample_chunks( pose::Pose & pose,
																 utility::vector1< SWA_Move > & swa_moves );

	void
	get_potential_resample_chunks( pose::Pose & pose,
																 utility::vector1< SWA_Move > & swa_moves );

	void
	remove_from_consideration_first_multi_residue_chunk( utility::vector1< SWA_Move > & swa_moves,
																											 bool remove_even_if_not_singlet );

	void
	get_potential_terminal_chunks( pose::Pose & pose,
																 utility::vector1< SWA_Move > & swa_moves,
																 AddOrDeleteChoice const & choice );

	void
	get_potential_add_chunks( pose::Pose & pose,
														utility::vector1< SWA_Move > & swa_moves );

	void
	get_random_chunk_at_chain_terminus( pose::Pose & pose,
																			SWA_Move & swa_move,
																			bool const disallow_delete,
																			bool const disallow_resample,
																			utility::vector1< Size > const & sample_res );

	void
	get_random_chunk_at_chain_terminus( pose::Pose & pose,
																		  SWA_Move & swa_move,
																			bool const disallow_delete,
																			bool const disallow_resample = true );

	void
	filter_by_sample_res(utility::vector1< SWA_Move > & swa_moves,
											 utility::vector1< Size > const & sample_res,
											 utility::vector1< Size > const & res_list );


} // monte_carlo
} // swa
} // protocols

#endif
