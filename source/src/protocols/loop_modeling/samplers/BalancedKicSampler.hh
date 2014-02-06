// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_samplers_BalancedKicSampler_HH
#define INCLUDED_protocols_loop_modeling_samplers_BalancedKicSampler_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/samplers/BalancedKicSampler.fwd.hh>
#include <protocols/loop_modeling/loggers/Logger.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/kinematic_closure/BalancedKicMover.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace samplers {

using protocols::kinematic_closure::perturbers::PerturberOP;
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;

class BalancedKicSampler : public LoopMoverTask {

public:
	BalancedKicSampler();
	BalancedKicSampler(loggers::LoggerOP logger);
	string get_name() const { return "BalancedKicSampler"; }

public:
	bool apply(Pose & pose, Loop const & loop, ScoreFunctionCOP score_function);

public:
	void add_perturber(PerturberOP perturber);
	void set_pivot_picker(PivotPickerOP picker);
	void log_filters(loggers::LoggerOP logger);

private:
	protocols::kinematic_closure::BalancedKicMover mover_;
};

}
}
}

#endif

