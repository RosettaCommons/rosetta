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
	get_random_residue_at_chain_terminus( pose::Pose & pose,
																				Size & residue_at_chain_terminus,
																				MovingResidueCase & moving_residue_case,
																				AddOrDeleteChoice & add_or_delete_choice,
																				bool const disallow_delete );

} // monte_carlo
} // swa
} // protocols

#endif
