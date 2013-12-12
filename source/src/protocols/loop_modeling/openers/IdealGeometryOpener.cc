// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/openers/IdealGeometryOpener.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/util.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

namespace protocols {
namespace loop_modeling {
namespace openers {

void IdealGeometryOpener::apply(Pose & pose, Loop const & loop) {
	for (Size index = loop.start(); index <= loop.stop(); index++) {
		core::conformation::idealize_position(index, pose.conformation());
	}
}

}
}
}

