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

/// @file   test/utility/signals/TokenHub.cxxtest.hh
/// @brief  test for utility::signals::TokenHub
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <utility/signals/TokenHub.hh>

// C++ headers
#include <vector>


// here must define a different namespace than in SignalHub.cxxtest.hh/BufferedSignalHub.cxxtest.hh
// otherwise we run into strange crosstalk situation and tests begin to fail
namespace sig4 {

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

} // namespace sig4


// --------------- Test Class --------------- //

class TokenHubTests : public CxxTest::TestSuite {


public: // setup


  // shared initialization
  void setUp() {
  }


  // shared finalization
  void tearDown() {
  }


  // --------------- Test Cases --------------- //

	/// @brief test TokenHub
	void test_TokenHub() {
		using namespace sig4;
		using namespace utility::signals;

		TokenHub< void, Event > hub1, hub2;

		Observer o1, o2, o3;

		Link l1 = hub1.connect_tokenize( &Observer::on_update, &o1 );
		Link l2 = hub1.connect_tokenize( &Observer::on_update, &o2 );
		Link l3 = hub1.connect_tokenize( &Observer::on_update, &o3 );

		TS_ASSERT_EQUALS( 3, hub1.n_tokens() );

		hub2.receive_tokens_from( hub1 );
		TS_ASSERT_EQUALS( 0, hub1.n_tokens() );
		TS_ASSERT_EQUALS( 3, hub2.n_tokens() );

		Event e( 100 );
		hub2( e );

		TS_ASSERT_EQUALS( 100, o1.recent_event.id );
		TS_ASSERT_EQUALS( 100, o2.recent_event.id );
		TS_ASSERT_EQUALS( 100, o3.recent_event.id );

		l2.invalidate();

		{ // begin scope
			TokenHub< void, Event > hub3;
			hub3.receive_tokens_from( hub2 );
			TS_ASSERT_EQUALS( 0, hub2.n_tokens() );
			TS_ASSERT_EQUALS( 2, hub3.n_tokens() );

			Event e2( 200 );
			hub3( e2 );

			TS_ASSERT_EQUALS( 200, o1.recent_event.id );
			TS_ASSERT_EQUALS( 100, o2.recent_event.id );
			TS_ASSERT_EQUALS( 200, o3.recent_event.id );

			hub3.disconnect( &Observer::on_update, &o3 );

			TS_ASSERT( !l3.valid() );

			Event e3( 300 );
			hub3( e3 );

			TS_ASSERT_EQUALS( 300, o1.recent_event.id );
			TS_ASSERT_EQUALS( 100, o2.recent_event.id );
			TS_ASSERT_EQUALS( 200, o3.recent_event.id );
		} // end scope

		TS_ASSERT( !l1.valid() );
	}


}; // class BufferedSignalHubTests

