// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/loggers/ScoreVsRmsd.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/moves/MonteCarlo.hh>

// Utility headers
#include <boost/foreach.hpp>

// C++ includes
#include <iostream>
#include <iomanip>
#include <algorithm>

#define foreach BOOST_FOREACH

namespace protocols {
namespace loop_modeling {
namespace loggers {

using namespace std;
using protocols::moves::MonteCarlo;

ScoreVsRmsd::ScoreVsRmsd(Pose const & native, Loop const & loop) {
	native_ = native;
	loops_ = Loops();
	loops_.add_loop(loop);
}

void ScoreVsRmsd::log_monte_carlo_(MonteCarlo const & monte_carlo) {
	using protocols::loops::loop_rmsd;

	if (monte_carlo.mc_accepted()) {
		Pose const & pose = monte_carlo.last_accepted_pose();
		Real score = monte_carlo.last_accepted_score();
		Real rmsd = loop_rmsd(native_, pose, loops_);

		map_[score] = rmsd;
	}
}

void ScoreVsRmsd::log_ending_(Pose const &) {
	ofstream file("score_vs_rmsd.txt");

	file << fixed;

	foreach (ScoreToRmsdMap::value_type const & pair, map_) {
		file << setw(8)  << setprecision(4) << pair.second << " ";	// RMSD
		file << setw(10) << setprecision(3) << pair.first << endl;	// Score
	}

	file.close();
}

}
}
}
