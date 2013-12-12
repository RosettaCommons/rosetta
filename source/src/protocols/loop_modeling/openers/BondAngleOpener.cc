// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/openers/types.hh>
#include <protocols/loop_modeling/openers/BondAngleOpener.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

// Numeric headers
#include <numeric/random/random.hh>

namespace protocols {
namespace loop_modeling {
namespace openers {

void BondAngleOpener::apply(Pose & pose, Loop const & loop) {
	using core::id::AtomID;
	using core::conformation::Conformation;
	using numeric::random::gaussian;

	Conformation & conformation = pose.conformation();
	AtomID ids[4];

	for (Size index = loop.start(); index <= loop.stop(); index++) {
		// Extracted from Top8000
		Real bond_angle = 1.9415 + 0.0398 * gaussian();

		ids[0] = AtomID(1, index);
		ids[1] = AtomID(2, index);
		ids[2] = AtomID(3, index);
		
		conformation.set_bond_angle(ids[0], ids[1], ids[2], bond_angle);
	}
}

}
}
}

