// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Kale Kundert (kale.kundert@ucsf.edu)

// Headers {{{1
#include <protocols/canonical_sampling/ProgressBarObserver.hh>
#include <protocols/canonical_sampling/ProgressBarObserverCreator.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// }}}1

namespace protocols {
namespace canonical_sampling {

// Global names {{{1
using namespace std;
using protocols::moves::MoverOP;

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.ProgressBar" );
// }}}1

string ProgressBarObserverCreator::keyname() const { // {{{1
	return ProgressBarObserverCreator::mover_name();
}

MoverOP ProgressBarObserverCreator::create_mover() const { // {{{1
	return MoverOP( new ProgressBarObserver );
}

string ProgressBarObserverCreator::mover_name() { // {{{1
	return "ProgressBarObserver";
}
// }}}1

ProgressBarObserver::ProgressBarObserver() // {{{1
: ThermodynamicObserver(), progress_(0) {}

ProgressBarObserver::ProgressBarObserver( // {{{1
	ProgressBarObserver const & other)
: ThermodynamicObserver(other), progress_(other.progress_) {}

ProgressBarObserver::~ProgressBarObserver() {} // {{{1

MoverOP ProgressBarObserver::clone() const { // {{{1
	return MoverOP( new ProgressBarObserver( *this ) );
}

void ProgressBarObserver::observe_after_metropolis( // {{{1
	MetropolisHastingsMover const & mover) {

	cout << "\r[" << ++progress_ << "/" << mover.ntrials() << "]" << flush;
	if ( progress_ == mover.ntrials() ) cout << endl;
}
// }}}1

} // namespace moves
} // namespace protocols

