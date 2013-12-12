// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/samplers/Opener.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocols headers
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loop.hh>

// Utility headers
#include <boost/utility.hpp>

namespace protocols {
namespace loop_modeling {
namespace samplers {

using protocols::moves::MoverOP;
using protocols::loops::Loop;

MoverOP Opener::setup(Loop const & loop) {
	loop_ = loop;
	return this;
}

}
}
}

