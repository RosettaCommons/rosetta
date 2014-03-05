// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/loggers/Logger.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.hh>

#include <utility/string_util.hh>

// C++ headers
#include <sstream>
#include <iomanip>
#include <cmath>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using namespace std;
using protocols::moves::MonteCarlo;
using protocols::moves::MonteCarloCOP;
using utility::vector1;

void Logger::log_beginning(Pose const & pose, Size I, Size J, Size K) {
	i_ = 1; j_ = 1; k_ = 1;
	I_ = I; J_ = J; K_ = K;
	log_beginning_(pose);
}

void Logger::log_beginning(Pose const & pose, vector1<Size> iterations) {
	log_beginning(pose, iterations[1], iterations[2], iterations[3]);
}

void Logger::log_iteration(Pose const & pose, Size i, Size j, Size k) {
	i_ = i; j_ = j; k_ = k;
	log_iteration_(pose);
}

void Logger::log_monte_carlo(MonteCarlo const & monte_carlo) {
	log_monte_carlo_(monte_carlo);
}

void Logger::log_monte_carlo(MonteCarloCOP monte_carlo) {
	log_monte_carlo_(*monte_carlo);
}

void Logger::log_ending(Pose const & pose)  {
	log_ending_(pose);
}

/// @details All of the index variables (i, j, k) count from 1.  This function
/// also counts from 1, in the sense that it will return 1 when i = j = k = 1.
Size Logger::get_iteration_as_int() const {
	return (i_ - 1) * J_ * K_ + (j_ - 1) * K_ + (k_ - 1) + 1;
}

Size Logger::get_max_iteration_as_int() const {
	return I_ * J_ * K_;
}

string Logger::get_iteration_as_string() const {
	Size iteration = get_iteration_as_int();
	Size max_iteration = get_max_iteration_as_int();
	Size max_digits = get_num_digits(max_iteration);

	stringstream stream;
	stream << setw(max_digits) << setfill('0') << iteration;
	return stream.str();
}

string Logger::get_max_iteration_as_string() const {
	stringstream stream;
	stream << get_max_iteration_as_int();
	return stream.str();
}

Size Logger::get_num_digits(Size value) const {
	return utility::get_num_digits( value );
}

}
}
}

