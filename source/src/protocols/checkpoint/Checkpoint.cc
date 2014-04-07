// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
#ifdef BOINC
#include <protocols/boinc/boinc.hh>
#endif // BOINC

// Unit headers
#include <protocols/checkpoint/Checkpoint.hh>



namespace protocols {
namespace checkpoint {

	void checkpoint_with_interval( const int interval_in ) {
		Timer timer = Timer::instance();
		timer.set_interval( interval_in );
	}

	Timer& Timer::instance(void)
	{
		static Timer theTimer;
		return theTimer;
	}

	void Timer::set_interval( const int interval_in ) {
		if (interval_in < 1) return;
		interval_ = interval_in;
	}

	bool Timer::is_on() { return is_on_; }

	bool Timer::time_to_checkpoint() {
#ifdef BOINC
		// honor user's disk write interval preference
		// note: boinc_time_to_checkpoint sets critical section to true if
		// it's time to checkpoint to prevent the client from stopping the app
		// when checkpointing. uses boinc api calls

		if (!boinc_time_to_checkpoint() && !boinc_is_standalone()) return false;
		boinc_end_critical_section(); // boinc_time_to_checkpoint sets is, and we dont need it yet.
#endif
		double time_diff(0.0);
		time_t curr_time;
		time(&curr_time);
		time_diff = difftime( curr_time, time_ );
		if ( time_diff > interval_ ) return true;

		return false;
	}

	void Timer::reset() {
		time(&time_);
#ifdef BOINC
		protocols::boinc::Boinc::update_pct_complete();
		protocols::boinc::Boinc::checkpoint_completed();
#endif
	}

	bool Timer::is_on_ = false;
	int Timer::interval_ = 600;
	time_t Timer::time_ = time(0);


} // checkpoint

} // protocols


