// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Headers {{{1
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/LoopProtocolCreator.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.hh>
#include <protocols/loop_modeling/utilities/AcceptanceCheck.hh>

// Protocols headers
#include <protocols/filters/Filter.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector1.hh>

// C++ headers
#include <iostream>
#include <cmath>

#define foreach BOOST_FOREACH

// Namespaces {{{1
using namespace std;

using core::Real;
using core::Size;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::moves::MoverOP;
using protocols::moves::MonteCarloOP;
using utility::tools::make_vector1;

using utility::tag::TagCOP;
using basic::datacache::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
// }}}1

namespace protocols {
namespace loop_modeling {

MoverOP LoopProtocolCreator::create_mover() const { // {{{1
	return MoverOP( new LoopProtocol );
}

string LoopProtocolCreator::keyname() const { // {{{1
	return "LoopProtocol";
}
// }}}1

LoopProtocol::LoopProtocol() { // {{{1
	protocol_ = utility::pointer::static_pointer_cast< utilities::LoopMoverGroup >( register_nested_loop_mover( LoopMoverOP( new utilities::LoopMoverGroup ) ) );
	movers_ = protocol_->add_mover_group();
	refiners_ = protocol_->add_mover_group();
	monte_carlo_ = NULL;
	sfxn_cycles_ = 1;
	temp_cycles_ = 1;
	mover_cycles_ = 1;
	ramp_sfxn_rep_ = false;
	ramp_sfxn_rama_ = false;
	ramp_temp_ = true;
	scale_temp_cycles_ = false;
	initial_temp_ = 1.5;
	final_temp_ = 0.5;
}

LoopProtocol::~LoopProtocol() {} // {{{1

void LoopProtocol::parse_my_tag( // {{{1
		TagCOP tag,
		DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const & pose) {

	sfxn_cycles_ = tag->getOption<Size>("sfxn_cycles", sfxn_cycles_);
	mover_cycles_ = tag->getOption<Size>("mover_cycles", mover_cycles_);
	ramp_sfxn_rama_ = tag->getOption<bool>("ramp_rama", ramp_sfxn_rama_);
	ramp_sfxn_rep_ = tag->getOption<bool>("ramp_rep", ramp_sfxn_rep_);
	ramp_temp_ = tag->getOption<bool>("ramp_temp", ramp_temp_);
	initial_temp_ = tag->getOption<Real>("initial_temp", initial_temp_);
	final_temp_ = tag->getOption<Real>("final_temp", final_temp_);

	// The temp_cycles tag is complicated because it's a number optionally 
	// followed by an 'x'.  If the 'x' is present, it means the number of cycles 
	// should be the specified number multiplied by the loop length.  Otherwise 
	// the number of cycles will just be the specified number.

	if (tag->hasOption("temp_cycles")) {
		using namespace boost::xpressive;

		string cycles = tag->getOption<string>("temp_cycles");
		sregex regex = (s1 = +digit) >> (s2 = !as_xpr('x'));
		smatch match;

		if (regex_match(cycles, match, regex)) {
			temp_cycles_ = boost::lexical_cast<Size>(match[1]);
			scale_temp_cycles_ = (match[2] == "x");
		} else {
			stringstream message;
			message << "getOption: key = temp_cycles stream extraction failed! ";
			message << "Tried to parse '" << cycles << "'\n";
			throw utility::excn::EXCN_Msg_Exception(message.str());
		}
	}

	// By default, the standard set of refiners will be applied after any 
	// manually specified movers.  This makes it easy to experiment with 
	// different movers without messing up the core protocol.  However, if you 
	// want to specify your own refinement algorithm, you need to disable this 
	// behavior via the auto_refine flag.
	
	if (! tag->getOption<bool>("auto_refine", true)) {
		clear_refiners();
	}

	// Add loop movers to this protocol by parsing the subtags.

	foreach (TagCOP subtag, tag->getTags()) {
		LoopMoverOP loop_mover = utilities::loop_mover_from_tag(subtag, data, filters, movers, pose);
		add_mover(loop_mover);
	}
}
// }}}1

bool LoopProtocol::do_apply(Pose & pose) { // {{{1
	start_protocol(pose);

	Size temp_cycles = scale_temp_cycles_ ? 
		temp_cycles_ * get_loops().loop_size() :
		temp_cycles_;

	for (Size i = 1; i <= sfxn_cycles_; i++) {
		monte_carlo_->recover_low(pose);
		//ramp_score_function(i);//commented out due to division by zero error when running the integration test; needs debugging

		for (Size j = 1; j <= temp_cycles; j++) {
			ramp_temperature(j);

			for (Size k = 1; k <= mover_cycles_; k++) {
				attempt_loop_move(pose, i, j, k);
			}
		}
		// The legacy code repacks and minimizes here.  This doesn't really 
		// translate to the new framework, because the "kic" loop mover isn't 
		// distinct from the "repack" and "minimize" loop movers.
	}

	pose = monte_carlo_->lowest_score_pose();
	finish_protocol(pose);

	return true;
}

void LoopProtocol::start_protocol(Pose & pose) { // {{{1
	if (movers_->empty() && refiners_->empty()) {
		utility_exit_with_message("Refusing to run LoopProtocol without any movers or refiners.");
	}

	// Setup the Monte Carlo simulation.

	protocols::loops::add_cutpoint_variants(pose);
	protocols::loops::loop_mover::loops_set_chainbreak_weight(
			get_score_function(), 1);

	monte_carlo_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo(
			pose, *get_score_function(), initial_temp_) );

	// Setup the loggers.

	BOOST_FOREACH (loggers::LoggerOP logger, loggers_) {
		logger->log_beginning(
				pose, make_vector1(sfxn_cycles_, temp_cycles_, mover_cycles_));
	}
}

void LoopProtocol::ramp_score_function(Size iteration) { // {{{1
	using core::scoring::fa_rep;
	using core::scoring::rama;

	Real ramp_factor = 1 / (iteration - sfxn_cycles_ + 1);//this line can lead to a division by zero exception

	if (ramp_sfxn_rep_) {
		Real repulsive_weight = get_score_function()->get_weight(fa_rep);
		get_score_function()->set_weight(fa_rep, repulsive_weight * ramp_factor);
	}
	if (ramp_sfxn_rama_) {
		Real rama_weight = get_score_function()->get_weight(rama);
		get_score_function()->set_weight(rama, rama_weight * ramp_factor);
	}
}

void LoopProtocol::ramp_temperature(Size /*iteration*/) { // {{{1
	if (ramp_temp_) {
		Real temperature = monte_carlo_->temperature();
		Real ramp_factor = std::pow(
				final_temp_ / initial_temp_, 1.0 / (sfxn_cycles_ * temp_cycles_));
		monte_carlo_->set_temperature(temperature * ramp_factor);
	}
}

void LoopProtocol::attempt_loop_move( // {{{1
		Pose & pose, Size i, Size j, Size k) {

	BOOST_FOREACH(loggers::LoggerOP logger, loggers_) {
		 logger->log_iteration(pose, i, j, k);
	}

	protocol_->apply(pose);

	if (protocol_->was_successful()) {
		monte_carlo_->boltzmann(pose);
	} else {
		monte_carlo_->set_last_accepted_pose(pose);
	}

	BOOST_FOREACH(loggers::LoggerOP logger, loggers_) {
		 logger->log_monte_carlo(monte_carlo_);
	}
}

void LoopProtocol::finish_protocol(Pose & pose) { // {{{1

	// Clean up the loggers.
	BOOST_FOREACH (loggers::LoggerOP logger, loggers_) {
		 logger->log_ending(pose);
	}
}
// }}}1

void LoopProtocol::add_mover(LoopMoverOP mover) { // {{{1
	movers_->add_mover(mover);
}

void LoopProtocol::add_refiner(LoopMoverOP refiner) { // {{{1
	refiners_->add_mover(refiner);
}

void LoopProtocol::add_filter(protocols::filters::FilterOP filter) { // {{{1
	movers_->add_filter(filter);
}

void LoopProtocol::add_acceptance_check(string name) { // {{{1
	movers_->add_mover(LoopMoverOP( new utilities::AcceptanceCheck(monte_carlo_, name) ));
}


void LoopProtocol::add_logger(loggers::LoggerOP logger) { // {{{1
	loggers_.push_back(logger);
}

void LoopProtocol::clear_movers() { // {{{1
	movers_->clear();
}

void LoopProtocol::clear_refiners() { // {{{1
	refiners_->clear();
}

void LoopProtocol::clear_movers_and_refiners() { // {{{1
	clear_movers();
	clear_refiners();
}

void LoopProtocol::mark_as_default() { // {{{1
	movers_->mark_as_default();
	refiners_->mark_as_default();
}
// }}}1

void LoopProtocol::set_sfxn_cycles(Size x) { // {{{1
	sfxn_cycles_ = x;
}

void LoopProtocol::set_temp_cycles(Size x, bool times_loop_length) { // {{{1
	temp_cycles_ = x;
	scale_temp_cycles_ = times_loop_length;
}

void LoopProtocol::set_mover_cycles(Size x) { // {{{1
	mover_cycles_ = x;
}

void LoopProtocol::set_temperature_schedule(Real start, Real stop) { // {{{1
	initial_temp_ = start;
	final_temp_ = stop;
}

void LoopProtocol::set_temperature_ramping(bool value) { // {{{1
	ramp_temp_ = value;
}

void LoopProtocol::set_repulsive_term_ramping(bool value) { // {{{1
	ramp_sfxn_rep_ = value;
}

void LoopProtocol::set_rama_term_ramping(bool value) { // {{{1
	ramp_sfxn_rama_ = value;
}
// }}}1

}
}
