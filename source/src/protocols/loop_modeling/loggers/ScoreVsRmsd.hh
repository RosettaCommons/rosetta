// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_loggers_ScoreVsRmsd_HH
#define INCLUDED_protocols_loop_modeling_loggers_ScoreVsRmsd_HH

// Unit headers
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/loggers/ScoreVsRmsd.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

// C++ headers
#include <fstream>
#include <string>
#include <vector>
#include <map>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using namespace std;
using core::Real;
using core::pose::Pose;
using protocols::moves::MonteCarlo;
using protocols::loops::Loops;
using protocols::loops::Loop;
typedef map<Real, Real> ScoreToRmsdMap;

class ScoreVsRmsd : public Logger {

public:
	ScoreVsRmsd(Pose const & native, Loop const & loop);

public:
	void log_iteration_(Pose const & pose);
	void log_monte_carlo_(MonteCarlo const & monte_carlo);
	void log_ending_(Pose const & pose);

private:
	Pose native_;
	Loops loops_;
	ScoreToRmsdMap map_;
	ofstream file_;

};

}
}
}

#endif

