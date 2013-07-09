// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_SWA_MonteCarloUtil.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_monte_carlo_RNA_SWA_MonteCarloUtil_hh
#define INCLUDED_protocols_swa_monte_carlo_RNA_SWA_MonteCarloUtil_hh

#include <core/pose/Pose.fwd.hh>
#include <protocols/swa/monte_carlo/types.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace monte_carlo {

	MovingResidueCase
	get_moving_residue_case( pose::Pose const & pose, Size const i );

	void
	get_potential_delete_residues( pose::Pose & pose,
																 utility::vector1< Size > & possible_res,
																 utility::vector1< MovingResidueCase > & moving_residue_cases,
																 utility::vector1< AddOrDeleteChoice > & add_or_delete_choices );

	void
	get_potential_resample_residues( pose::Pose & pose,
																	 utility::vector1< Size > & possible_res );

	void
	get_potential_resample_residues( pose::Pose & pose,
																	 utility::vector1< Size > & possible_res,
																	 utility::vector1< MovingResidueCase > & moving_residue_cases,
																	 utility::vector1< AddOrDeleteChoice > & add_or_delete_choices );

	void
	get_potential_terminal_residues( pose::Pose & pose,
																	 utility::vector1< Size > & possible_res,
																	 utility::vector1< MovingResidueCase > & moving_residue_cases,
																	 utility::vector1< AddOrDeleteChoice > & add_or_delete_choices,
																	 AddOrDeleteChoice const & choice );

	void
	get_potential_add_residues( pose::Pose & pose,
															utility::vector1< Size > & possible_res,
															utility::vector1< MovingResidueCase > & moving_residue_cases,
															utility::vector1< AddOrDeleteChoice > & add_or_delete_choices );

	void
	get_random_residue_at_chain_terminus( pose::Pose & pose,
																				Size & residue_at_chain_terminus,
																				MovingResidueCase & moving_residue_case,
																				AddOrDeleteChoice & add_or_delete_choice,
																				bool const disallow_delete,
																				bool const disallow_resample = true );

} // monte_carlo
} // swa
} // protocols

#endif
