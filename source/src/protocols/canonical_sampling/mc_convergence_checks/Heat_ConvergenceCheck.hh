// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_Heat_ConvergenceCheck_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_Heat_ConvergenceCheck_hh


// type headers
#include <core/types.hh>

// unit headers
#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <utility/excn/Exceptions.hh>

// package headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

class EXCN_Heat_Converged : public moves::EXCN_Converged {};

class Heat_ConvergenceCheck : public moves::MonteCarloExceptionConverge {
	virtual bool operator() ( const core::pose::Pose&, moves::MonteCarlo const& mc, bool /*reject*/ ) {
		if ( mc.last_accept() >= mc.heat_after_cycles() - mc.check_frequency() ) {
			throw canonical_sampling::mc_convergence_checks::EXCN_Heat_Converged();
		}
		return true;
	}
private:
};

} // mc_convergence_check
} // moves
} // rosetta

#endif
