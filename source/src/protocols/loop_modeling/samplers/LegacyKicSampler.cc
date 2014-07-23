// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/samplers/LegacyKicSampler.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>

// C++ headers
#include <iostream>

using namespace std;

namespace protocols {
namespace loop_modeling {
namespace samplers {

using core::Size;
using core::pose::Pose;
using protocols::loops::loop_closure::kinematic_closure::KinematicMover;
using protocols::loops::loop_closure::kinematic_closure::TorsionSamplingKinematicPerturber;
using protocols::loops::loop_closure::kinematic_closure::TorsionSamplingKinematicPerturberOP;

LegacyKicSampler::LegacyKicSampler() { // {{{1

	mover_ = new KinematicMover();
	TorsionSamplingKinematicPerturberOP perturber =
		new TorsionSamplingKinematicPerturber(mover_.get());

	//perturber->set_max_sample_iterations(1);
	mover_->set_perturber(perturber);
}

bool LegacyKicSampler::do_apply(Pose & pose, Loop const & loop) { // {{{1
	Size pivot_1 = loop.start();
	Size pivot_3 = loop.stop();
	Size pivot_2 = pivot_1 + (pivot_3 - pivot_1) / 2;

	mover_->set_pivots(pivot_1, pivot_2, pivot_3);
	mover_->apply(pose);

	return true;
}
// }}}1

}
}
}

