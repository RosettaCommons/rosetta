// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/MultiTempTrialCounter.cc
/// @brief MultiTempTrialCounter methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/MultiTempTrialCounter.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>

// Core Headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>

// External headers
#include <boost/foreach.hpp>

static basic::Tracer tr("protocols.canonical_sampling.MultiTempTrialCounter");

namespace protocols {
namespace canonical_sampling {

using namespace std;
using namespace core;
using namespace ObjexxFCL;
using protocols::moves::MoverOP;
using protocols::moves::TrialCounter;
using utility::vector1;

MultiTempTrialCounter::MultiTempTrialCounter(TemperatureControllerCOP tc) // {{{1
	: TrialCounter(), temp_controller_(tc) {
	reset();
}

MultiTempTrialCounter::~MultiTempTrialCounter() {} // {{{1
// }}}1

void MultiTempTrialCounter::reset() { // {{{1
	counters_.clear();
	counters_.resize(temp_controller_->n_temp_levels());
}

void MultiTempTrialCounter::reset_temp_controller( // {{{1
		TemperatureControllerCOP temp_controller) {
	temp_controller_ = temp_controller;
	reset();
}

void MultiTempTrialCounter::count_trial(string const & tag) { // {{{1
	get_current_counter().count_trial(tag);
}

void MultiTempTrialCounter::count_accepted(string const & tag) { // {{{1
	get_current_counter().count_accepted(tag);
}

void MultiTempTrialCounter::count_energy_drop(string const & tag, Real drop) { // {{{1
	get_current_counter().count_energy_drop(tag, drop);
}

void MultiTempTrialCounter::collect() { // {{{1
#ifdef USEMPI

	// The idea is that I want to be able to read all of the move statistics for 
	// the entire simulation out of a single trial counter object.  One reason is 
	// that I want to be able to do single-threaded output even if my protocol is 
	// multi-threaded.  The problem is that the trial counter from each thread 
	// collects its own statistics, and so no complete records exist.  The 
	// solution is to have this collect step.
	//
	// I'm envisioning that this method will begin with MPI_Gather so the root 
	// node can collect statistics from every thread.  Then MPI_Scatter will be 
	// called to broadcast the statistics to every thread.  By the end of the 
	// method, every counter will contain a complete record of every move.
	
	utility_exit_with_message("Not Implemented");
#endif
}

TrialCounter const & MultiTempTrialCounter::temp_level(Size level) const { // {{{1
	assert(temp_controller_->n_temp_levels() == counters_.size());
	return counters_[level];
}

Size MultiTempTrialCounter::num_temp_levels() const { // {{{1
	return counters_.size();
}

Size MultiTempTrialCounter::total_trials() const { // {{{1
	Size total_trials = 0;
	BOOST_FOREACH(TrialCounter const & counter, counters_) {
		total_trials += counter.total_trials();
	}
	return total_trials;
}

Size MultiTempTrialCounter::trial(std::string const & tag) const { // {{{1
	Size num_attempted = 0;
	BOOST_FOREACH(TrialCounter const & counter, counters_) {
		num_attempted += counter.trial(tag);
	}
	return num_attempted;
}

Size MultiTempTrialCounter::accepted(std::string const & tag) const { // {{{1
	Size num_accepted = 0;
	BOOST_FOREACH(TrialCounter const & counter, counters_) {
		num_accepted += counter.accepted(tag);
	}
	return num_accepted;
}

Real MultiTempTrialCounter::energy_drop(std::string const & tag) const { // {{{1
	Size total_drop = 0;
	BOOST_FOREACH(TrialCounter const & counter, counters_) {
		total_drop += counter.energy_drop(tag);
	}
	return total_drop;
}

vector1<string> const MultiTempTrialCounter::tags() const { // {{{1
	vector1<string> tags, subtags;
	BOOST_FOREACH(TrialCounter const & counter, counters_) {
		subtags = counter.tags();
		tags.insert(tags.end(), subtags.begin(), subtags.end());
	}
	return tags;
}

void MultiTempTrialCounter::show( // {{{1
		ostream& os, string, bool ) const
{
	write_to_stream(os, "");
}

void MultiTempTrialCounter::write_to_file( // {{{1
		string const & file, string const & tag) const {

	utility::io::ozstream out(file, std::ios::app);
	write_to_stream(out, tag);
}
// }}}1

TrialCounter & MultiTempTrialCounter::get_current_counter() { // {{{1
	assert(temp_controller_);
	assert(temp_controller_->n_temp_levels() == counters_.size());
	return counters_[temp_controller_->temperature_level()];
}

void MultiTempTrialCounter::write_to_stream( // {{{1
		std::ostream& os, std::string const& tag) const {

	assert(temp_controller_);
	for (Size i=1; i<=counters_.size(); ++i) {
		std::ostringstream line_header;
		line_header << "level " << format::I(4, i) << " temperature " << format::F(4, 2, temp_controller_->temperature(i));
		counters_[ i ].show(os, line_header.str(), false);
		os << " " << tag << std::endl;
	}
}
// }}}1

} //moves
} //protocols

