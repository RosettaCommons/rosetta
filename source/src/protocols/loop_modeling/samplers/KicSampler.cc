// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/samplers/KicSampler.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/solution_pickers/SolutionPicker.hh>

namespace protocols {
namespace loop_modeling {
namespace samplers {

using namespace std;

KicSampler::KicSampler() {}

KicSampler::KicSampler(loggers::LoggerOP logger) {
	log_filters(logger);
}

void KicSampler::setup(
		Pose & pose, Loop const & loop, ScoreFunctionOP) {

	mover_.setup(pose, loop);
}

bool KicSampler::apply(
		Pose & pose, Loop const & loop, ScoreFunctionCOP) {

	mover_.setup(pose, loop);
	mover_.apply(pose);
	return true;
}

void KicSampler::add_perturber(PerturberOP perturber) {
	mover_.add_perturber(perturber);
}

void KicSampler::set_pivot_picker(PivotPickerOP picker) {
	mover_.set_pivot_picker(picker);
}

void KicSampler::set_solution_picker(SolutionPickerOP picker) {
	mover_.set_solution_picker(picker);
}

void KicSampler::log_filters(loggers::LoggerOP logger) {
	mover_.log_filters(logger);
}

}
}
}

