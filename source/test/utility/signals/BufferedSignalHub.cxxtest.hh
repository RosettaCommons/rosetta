// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   test/utility/signals/BufferedSignalHub.cxxtest.hh
/// @brief  test for utility::signals::BufferedSignalHub
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <utility/signals/BufferedSignalHub.hh>

// C++ headers
#include <vector>


// here must define a different namespace than in SignalHub.cxxtest.hh
// otherwise we run into strange crosstalk situation and tests begin to fail
namespace sig2 {


// Event
struct Event {
	Event() {};
	Event( int i ) : id( i ) {}

	int id;
};


// Observer
struct Observer {

	Observer() : count( 0 ) {}

	void
	on_update( Event e ) {
		++count;
		recent_event = e;
	}

	unsigned long count;
	Event recent_event;

}; // class Observer


} // namespace sig2


// --------------- Test Class --------------- //

class BufferedSignalHubTests : public CxxTest::TestSuite {


public: // setup


  // shared initialization
  void setUp() {
  }


  // shared finalization
  void tearDown() {
  }


  // --------------- Test Cases --------------- //

	/// @brief test observer attachment, detachment, and notification
	void test_BufferedSignalHub_buffering() {
		using namespace sig2;

		utility::signals::BufferedSignalHub< void, Event > hub;
		Observer o;

		hub.connect( &Observer::on_update, &o );
		hub.buffer();

		TS_ASSERT( hub.blocked() );
		TS_ASSERT( hub.buffering() );

		// send 3 events
		hub( Event( 5 ) );
		hub( Event( 6 ) );
		hub( Event( 7 ) );

		TS_ASSERT_EQUALS( hub.buffer_size(), 3 );

		// unblock and release the buffer
		hub.unblock();

		TS_ASSERT_EQUALS( o.count, 3 );
		TS_ASSERT_EQUALS( o.recent_event.id, 7 );

		TS_ASSERT_EQUALS( hub.buffer_size(), 0 );
		TS_ASSERT( !hub.blocked() );
		TS_ASSERT( !hub.buffering() );
	}


}; // class BufferedSignalHubTests

