// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/samplers/BalancedKicSampler.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/kinematic_closure/BalancedKicMover.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace loop_modeling {
namespace samplers {

using namespace std;

BalancedKicSampler::BalancedKicSampler() {} 

BalancedKicSampler::BalancedKicSampler(loggers::LoggerOP logger) {
	log_filters(logger);
}

bool BalancedKicSampler::apply(
		Pose & pose, Loop const & loop, ScoreFunctionCOP score_function) {

	mover_.set_loop(loop);
	mover_.apply(pose);
	return true;
}

void BalancedKicSampler::add_perturber(PerturberOP perturber) {
	mover_.add_perturber(perturber);
}

void BalancedKicSampler::set_pivot_picker(PivotPickerOP picker) {
	mover_.set_pivot_picker(picker);
}

void BalancedKicSampler::log_filters(loggers::LoggerOP logger) {
	mover_.log_filters(logger);
}

}
}
}

