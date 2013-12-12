// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/openers/RamachandranOpener.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

namespace protocols {
namespace loop_modeling {
namespace openers {

void RamachandranOpener::apply(Pose & pose, Loop const & loop) {

	// This function has been simplified in many ways compared to the standard 
	// KIC implementation.  First, no move-map is respected.  Second, no omega 
	// sampling is done under any circumstances.

	using core::scoring::Ramachandran;
	using core::scoring::ScoringManager;
	
	Real phi, psi;
	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	for (Size index = loop.start(); index <= loop.stop(); index++) {
		rama.random_phipsi_from_rama(pose.aa(index), phi, psi);
		pose.set_phi(index, phi); pose.set_psi(index, psi);
	}
}

}
}
}

