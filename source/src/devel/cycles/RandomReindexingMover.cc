// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <devel/cycles/RandomReindexingMover.hh>
#include <devel/cycles/ReindexingMover.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <numeric/random/random.hh>

namespace devel {
namespace cycles {

using namespace core;

/// A random offset will be chosen by picking a uniform integer between 0 and 
/// the length of the given pose.  The pose should have been previously passed 
/// through the SetupMover, otherwise this method may trigger some assertion 
/// failures.

void RandomReindexingMover::apply(pose::Pose &pose) {
	Size cycle_length = pose.size() - 1;
	Size offset = numeric::random::random_range(0, cycle_length);
	ReindexingMover reindexer(offset);
	reindexer.apply(pose);
}

} // End 'cycles' namespace.
} // End 'devel' namespace.


