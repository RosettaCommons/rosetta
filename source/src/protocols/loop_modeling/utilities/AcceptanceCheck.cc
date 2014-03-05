// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/utilities/AcceptanceCheck.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

using protocols::moves::MonteCarloOP;

AcceptanceCheck::AcceptanceCheck(MonteCarloOP monte_carlo, string name)
	: monte_carlo_(monte_carlo), name_(name) {}

bool AcceptanceCheck::do_apply(Pose & pose) {
	return monte_carlo_->boltzmann(pose, name_);
}

} // namespace utilities
} // namespace kinematic_closure
} // namespace protocols

