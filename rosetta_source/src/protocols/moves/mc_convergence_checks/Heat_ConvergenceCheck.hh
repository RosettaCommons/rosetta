// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_moves_mc_convergence_checks_Heat_ConvergenceCheck_hh
#define INCLUDED_protocols_moves_mc_convergence_checks_Heat_ConvergenceCheck_hh


// type headers
#include <core/types.hh>

// unit headers
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <utility/excn/Exceptions.hh>

// package headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace moves {
namespace mc_convergence_checks {

class EXCN_Heat_Converged : public EXCN_Converged {};

class Heat_ConvergenceCheck : public ConvergenceCheck {
	virtual bool operator() ( const core::pose::Pose&, moves::MonteCarlo const& mc, bool /*reject*/ ) {
		if ( mc.last_accept() >= mc.heat_after_cycles() - mc.check_frequency() ) {
			throw moves::mc_convergence_checks::EXCN_Heat_Converged();
		}
		return true;
	}
private:
};

} // mc_convergence_check
} // moves
} // rosetta

#endif
