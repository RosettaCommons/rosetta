// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/loggers/AcceptanceRates.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.hh>

// C++ includes
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using namespace std;

void AcceptanceRates::log_task_(Pose const &, string name, bool successful) {
	vector<string>::const_iterator begin = task_names_.begin();
	vector<string>::const_iterator end = task_names_.end();
	uint index = find(begin, end, name) - begin;
	bool name_found = (index < task_names_.size());

	if (name_found == false) {
		task_names_.push_back(name);
		task_calls_.push_back(0);
		task_successes_.push_back(0);
	}

	task_calls_[index] += 1;
	task_successes_[index] += successful;
}

void AcceptanceRates::log_monte_carlo_(MonteCarlo const & monte_carlo) {
	if (all_tasks_succeeded()) {
		Pose const & pose = monte_carlo.last_accepted_pose();
		bool accepted = monte_carlo.mc_accepted() ? true : false;
		log_task_(pose, "MonteCarlo", accepted);
	}
}

void AcceptanceRates::log_ending_(Pose const &) {
	int name_chars = 0;
	int most_calls = 0;
	int most_successes = 0;

	for (uint i = 0; i < task_names_.size(); ++i) {
		name_chars = max(name_chars, (int) task_names_[i].size());
		most_calls = max(most_calls, task_calls_[i]);
		most_successes = max(most_successes, task_successes_[i]);
	}

	int call_digits = get_num_digits(most_calls);
	int success_digits = get_num_digits(most_successes);
	string underline(name_chars + call_digits + success_digits + 3, '=');
	
	cout << endl;
	cout << "Acceptance Rates" << endl;
	cout << underline << endl;

	for (uint i = 0; i < task_names_.size(); ++i) {
		string name = task_names_[i];

		cout << setw(name_chars + 2) << task_names_[i] + ": ";
		cout << setw(success_digits) << task_successes_[i] << "/";
		cout << setw(call_digits) << right << task_calls_[i] << endl;
	}
}

}
}
}

