// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/samplers/RamaOpener.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

// Numeric headers
#include <numeric/random/random.hh>

namespace protocols {
namespace loop_modeling {
namespace samplers {

using namespace std;
using core::Size;
using core::Real;
using core::pose::Pose;
using protocols::loops::LoopOP;

void RamaOpener::apply(Pose & pose) {

	// This function has been simplified in many ways compared to the standard 
	// KIC implementation.  First, no move-map is respected.  Second, no omega 
	// sampling is done under any circumstances.

	// Sample torsion angles:

	using core::scoring::Ramachandran;
	using core::scoring::ScoringManager;
	
	Real phi, psi;
	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	for (Size index = loop_.start(); index <= loop_.stop(); index++) {
			rama.random_phipsi_from_rama(pose.aa(index), phi, psi);
			pose.set_phi(index, phi); pose.set_psi(index, psi);
	}

	// Sample bond angles:

	/*
	using core::id::AtomID;
	using core::conformation::Conformation;
	using numeric::random::gaussian;

	Conformation & conformation = pose.conformation();
	AtomID ids[4];

	for (Size index = loop_.start(); index <= loop_.stop(); index++) {
		Real bond_angle = 111.24 + 2.28 * gaussian();

		ids[0] = AtomID(1, index);
		ids[1] = AtomID(2, index);
		ids[2] = AtomID(3, index);
		
		conformation.set_bond_angle(ids[0], ids[1], ids[2], bond_angle);
	}
	*/
}

}
}
}

