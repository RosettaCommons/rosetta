// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MoverStatus.hh
/// @brief  return status enum for Movers
/// @author Steven Lewis smlewi@gmail.com
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_protocols_moves_MoverStatus_hh
#define INCLUDED_protocols_moves_MoverStatus_hh

#include <string>

namespace protocols {
namespace moves {

/// @brief return status for movers - mover was successful, failed but can be retried, etc; used mostly by job dist.
/// @details Documention for individual codes from SML, and reflects the original intent behind them (more than current usage)
enum MoverStatus {
	/// @brief MS_SUCCESS: job succeeds.
	MS_SUCCESS = 0,
	/// @breif FAIL_RETRY: The job failed in a "science-y" way (rather than a "programmatic" way).  E.g. this particular
	/// docking run failed an RMSD filter, but there's no reason to believe it's not just trapped in a bad region of
	/// conformational space due to the stochasticity of Monte Carlo.
	FAIL_RETRY,
	/// @brief FAIL_DO_NOT_RETRY: This particular job should not be re-attempted, but there wasn't anything structurally
	/// wrong with it - it was a valid job, but its output is uninteresting, so toss it and move on.
	FAIL_DO_NOT_RETRY,
	/// @brief FAIL_BAD_INPUT: This job has something structurally wrong with it/its input and can never be completed.
	/// Other jobs of the same input (same InnerJob / stuff that differs only by its nstruct index) will by definition
	/// fail in the same fashion, so short-circuit those failures and pre-emptively cancel those jobs when possible.
	/// This probably means there is a user input problem.
	FAIL_BAD_INPUT,
	/// @brief FAIL: ??? Unknown usage (not obeyed by JD2).
	FAIL,


	//Alternative names of enums.
	MS_FAIL_RETRY = FAIL_RETRY,
	MS_FAIL_DO_NOT_RETRY = FAIL_DO_NOT_RETRY,
	MS_FAIL_BAD_INPUT = FAIL_BAD_INPUT,
	MS_FAIL = FAIL
};

MoverStatus mstype_from_name( std::string const & name );

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_MoverStatus_HH
