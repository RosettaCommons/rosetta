// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_loggers_Logger_HH
#define INCLUDED_protocols_loop_modeling_loggers_Logger_HH

// Unit headers
#include <protocols/loop_modeling/loggers/Logger.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <boost/utility.hpp>

// C++ headers
#include <string>

namespace protocols {
namespace loop_modeling {
namespace loggers {

using std::string;
using core::Size;
using core::pose::Pose;
using protocols::moves::MonteCarlo;
using protocols::moves::MonteCarloCOP;
using utility::pointer::ReferenceCount;
using boost::noncopyable;

class Logger : public ReferenceCount, protected noncopyable {

// Public Interface
public:

	void log_beginning(Pose const & pose, Size I, Size J, Size K);
	void log_beginning(Pose const & pose, utility::vector1<Size> iterations);
	void log_iteration(Pose const & pose, Size i, Size j, Size k);
	void log_task(Pose const & pose, string name, bool successful, bool required=true);
	void log_monte_carlo(MonteCarlo const & monte_carlo);
	void log_monte_carlo(MonteCarloCOP monte_carlo);
	void log_ending(Pose const & pose);

// Behavior Methods
protected:

	virtual void log_beginning_(Pose const & /*pose*/) {}
	virtual void log_iteration_(Pose const & /*pose*/) {}
	virtual void log_task_(Pose const & /*pose*/, string /*name*/, bool) {}
	virtual void log_monte_carlo_(MonteCarlo const & /*monte_carlo*/) {}
	virtual void log_ending_(Pose const & /*pose*/) {}

// Private Helpers
protected:

	Size get_iteration_as_int() const;
	Size get_max_iteration_as_int() const;
	string get_iteration_as_string() const;
	string get_max_iteration_as_string() const;
	Size get_num_digits(Size value) const;
	bool all_tasks_succeeded() const;

// Data Members
private:

	Size i_, j_, k_;		// Current iteration.
	Size I_, J_, K_;		// Maximum iterations.

	bool all_tasks_succeeded_;

};

}
}
}

#endif

