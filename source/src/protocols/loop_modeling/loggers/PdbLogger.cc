// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/loggers/PdbLogger.hh>

// Core includes
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.hh>

// C++ includes
#include <iostream>
#include <iomanip>
#include <cmath>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using namespace std;
using protocols::moves::MonteCarlo;

void PdbLogger::log_beginning_(Pose const & pose) {
	log_pose(pose, 0);
}

void PdbLogger::log_monte_carlo_(MonteCarlo const & monte_carlo) {
	log_pose(monte_carlo.last_accepted_pose());
}

void PdbLogger::log_pose(Pose const & pose) const {
	int iteration = get_iteration_as_int();
	log_pose(pose, iteration);
}

void PdbLogger::log_pose(Pose const & pose, Size iteration) const {
	int max_iteration = get_max_iteration_as_int();
	int max_digits = get_num_digits(max_iteration);

	stringstream stream;
	stream << "trajectory/" << setw(max_digits) << setfill('0') << iteration;
	stream << ".pdb";

	pose.dump_pdb(stream.str());
}

}
}
}
