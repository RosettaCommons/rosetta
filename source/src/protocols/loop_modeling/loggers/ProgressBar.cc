// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/loggers/ProgressBar.hh>

// Core includes
#include <core/pose/Pose.hh>

// Utility includes
#include <basic/Tracer.hh>

// C++ includes
#include <iostream>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using namespace std;
using core::pose::Pose;

static thread_local basic::Tracer TR( "protocols.loop_modeling.loggers.ProgressBar" );

ProgressBar::ProgressBar(string label) : label_(label) {}

void ProgressBar::log_iteration_(Pose const &) {
	int iteration = get_iteration_as_int();
	int max_iteration = get_max_iteration_as_int();

	TR << label_ << "[" << iteration << "/" << max_iteration << "]" << endl;
}

}
}
}

