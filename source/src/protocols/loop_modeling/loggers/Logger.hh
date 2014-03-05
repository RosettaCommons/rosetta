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
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/loggers/Logger.fwd.hh>

// Core headers
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

/// @brief Base class for logging information about the loop modeling run.
///
/// @details Subclasses should reimplement the protected method which all end 
/// with an underscore.  These methods are expected to do the actual logging, 
/// and will be called by the public methods of the same name.  In most cases, 
/// the public methods do a little bit of useful preprocessing before calling 
/// the corresponding protected methods.

class Logger :
	public utility::pointer::ReferenceCount, protected boost::noncopyable {

// Public Interface
public:

	/// @brief Record the beginning of a simulation.
	/// @details The arguments to this method should give the total number of 
	/// iterations that will be run in each loop of the simulation.
	void log_beginning(Pose const & pose, Size I, Size J, Size K);

	/// @brief Record the beginning of a simulation.
	/// @details The argument to this method should give the total number of 
	/// iterations that will be run in each loop of the simulation.
	void log_beginning(Pose const & pose, utility::vector1<Size> iterations);

	/// @brief Record the start of another iteration.
	/// @details The arguments to this method should identify the current 
	/// iteration.  Since the loop protocol utilizes three nested loops, three 
	/// index numbers are required.
	void log_iteration(Pose const & pose, Size i, Size j, Size k);

	/// @brief Record the result of a MonteCarlo acceptance check.
	void log_monte_carlo(protocols::moves::MonteCarlo const & monte_carlo);

	/// @brief Record the result of a MonteCarlo acceptance check.
	void log_monte_carlo(protocols::moves::MonteCarloCOP monte_carlo);

	/// @brief Record the end of a simulation.
	void log_ending(Pose const & pose);

// Behavior Methods
protected:

	/// @brief Reimplemented by subclasses to record the beginning of a 
	/// simulation.
	virtual void log_beginning_(Pose const &) {}

	/// @brief Reimplemented by subclasses to record the start of another 
	/// iteration.
	virtual void log_iteration_(Pose const &) {}

	/// @brief Reimplemented by subclasses to record the result of a MonteCarlo 
	/// acceptance check.
	virtual void log_monte_carlo_(protocols::moves::MonteCarlo const &) {}

	/// @brief Reimplemented by subclasses to record the end of a simulation.
	virtual void log_ending_(Pose const &) {}

// Private Helpers
protected:

	/// @brief Return the current iteration as an integer.
	Size get_iteration_as_int() const;

	/// @brief Return the number of iterations that will be run in the current 
	/// simulation as an integer.
	Size get_max_iteration_as_int() const;

	/// @brief Return the current iteration as a string.
	string get_iteration_as_string() const;

	/// @brief Return the number of iterations that will be run in the current 
	/// simulation as an string.
	string get_max_iteration_as_string() const;

	/// @brief Return the number of digits in the given number.
	/// @details This is useful for formatting operations.
	Size get_num_digits(Size value) const;

// Data Members
private:

	Size i_, j_, k_;		// Current iteration.
	Size I_, J_, K_;		// Maximum iterations.

};

}
}
}

#endif

