// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_loggers_AcceptanceRates_HH
#define INCLUDED_protocols_loop_modeling_loggers_AcceptanceRates_HH

// Unit headers
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/loggers/AcceptanceRates.fwd.hh>

// C++ headers
#include <vector>
#include <string>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using namespace std;

class AcceptanceRates : public Logger {

public:

	void log_task_(Pose const & pose, string name, bool successful);
	void log_monte_carlo_(MonteCarlo const & monte_carlo);
	void log_ending_(Pose const & pose);

private:

	vector<string> task_names_;
	vector<int> task_calls_;
	vector<int> task_successes_;

	bool all_tasks_succeeded_;

};

}
}
}

#endif

