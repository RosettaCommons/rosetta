// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_checkpoint_Checkpoint_hh
#define INCLUDED_protocols_checkpoint_Checkpoint_hh

#include <protocols/checkpoint/Checkpoint.fwd.hh>

// C++ Headers
#include <ctime>

namespace protocols {
namespace checkpoint {


void checkpoint_with_interval( const int interval_in );


///////////////////////////////////////////////////////////////////////////////
/// @brief: singleton checkpoint timer class
///
/// @details: Keeps track of when to checkpoint using a time interval.
///         Not thread safe.
///
///         This doesn't derive from SingletonBase because there's
///         currently no non-static data.
///
/// @author David K
///
///////////////////////////////////////////////////////////////////////////////
class Timer {


public:

	static Timer& instance(void);

	void set_interval( const int interval_in );
	static bool is_on(void);
	static bool time_to_checkpoint(void);
	static void reset(void);

private:
	Timer(void) {
		if ( !is_on_ ) time(&time_);
		is_on_ = true;
	}

	static bool is_on_; // do checkpointing?
	static int interval_; // checkpoint interval in seconds
	static time_t time_; // time of the last checkpoint

}; // Timer


} // checkpoint
} // protocols


#endif
