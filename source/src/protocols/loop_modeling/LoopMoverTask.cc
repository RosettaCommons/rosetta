// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/utilities/RepeatedTask.hh>
#include <protocols/loop_modeling/utilities/PeriodicTask.hh>

// Core headers
#include <core/types.hh>

namespace protocols {
namespace loop_modeling {

using core::Size;

LoopMoverTaskOP operator *(LoopMoverTaskOP task, Size iterations) {
	return new utilities::RepeatedTask(task, iterations);
}

LoopMoverTaskOP operator %(LoopMoverTaskOP task, Size period) {
	return new utilities::PeriodicTask(task, period);
}

}
}

