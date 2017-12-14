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

/// @file   test/utility/signals/SignalHub.cxxtest.hh
/// @brief  test for utility::signals::SignalHub
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <utility/signals/SignalHub.hh>
#include <utility/signals/Link.hh>

// C++ headers
#include <vector>


namespace sig {


// Event
struct Event {
	int type;
	int data;
};


// Subject
class Subject {


public:


	typedef utility::signals::SignalHub< void, sig::Event >::Size Size;
	typedef utility::signals::Link Link;


public:


	template< typename MemFn, typename Ptr >
	Link
	attach( MemFn fn, Ptr ptr ) {
		return hub_.connect( fn, ptr );
	}


	template< typename MemFn, typename Ptr >
	bool
	detach( MemFn fn, Ptr ptr ) {
		return hub_.disconnect( fn, ptr );
	}


	void
	notify( Event e ) {
		hub_( e );
	}


	Size
	n_obs() {
		return hub_.size();
	}


	void
	clear() {
		hub_.clear();
	}


private:
	utility::signals::SignalHub< void, sig::Event > hub_;


}; // class Subject


// Observer
class Observer {


public:
	Observer( int /*id*/ ) /*: id_( id )*/ {
	}


public:


	void
	on_update_1( Event e ) {
		event1_ = e;
	}


	void
	on_update_2( Event e ) {
		event2_ = e;
	}


	Event &
	event1() {
		return event1_;
	}


	Event &
	event2() {
		return event2_;
	}


private:

	Event event1_;
	Event event2_;

}; // class Observer


} // namespace sig


// --------------- Test Class --------------- //

class SignalHubTests : public CxxTest::TestSuite {


public: // setup


	// shared initialization
	void setUp() {
	}


	// shared finalization
	void tearDown() {
	}


	// --------------- Test Cases --------------- //

	/// @brief test observer attachment, detachment, and notification
	void test_SignalHub_attach_detach_notify() {
		using namespace sig;

		Subject s;
		std::vector< sig::Observer > v;

		// make 1000 Observers for testing
		for ( int i = 0; i < 1000; ++i ) {
			v.push_back( Observer( i ) );
		}

		// attach
		for ( std::vector< Observer >::size_type i = 0; i < 1000; ++i ) {
			s.attach( &Observer::on_update_1, &v[ i ] );
		}

		// remove all even observers
		for ( std::vector< Observer >::size_type i = 0; i < v.size(); ++i ) {
			if ( i % 2 == 0 ) {
				TS_ASSERT( s.detach( &Observer::on_update_1, &v[ i ] ) );
			}
		}

		// create event and notify
		Event e;
		e.type = 30;
		e.data = 1020;

		s.notify( e );
		for ( std::vector< Observer >::size_type i = 0; i < v.size(); ++i ) {
			if ( i % 2 == 1 ) { // only odd observers receive the event
				TS_ASSERT_EQUALS( 30, v[ i ].event1().type );
				TS_ASSERT_EQUALS( 1020, v[ i ].event1().data );
			}
		}

		// cleanup, remove all observers
		s.clear();
	}


	/// @brief trigger false attach/detach
	void test_SignalHub_bad_attach_detach() {
		using namespace sig;
		using utility::signals::Link;

		Subject s;

		// attached observer
		Observer o( 2000 );
		Link link = s.attach( &Observer::on_update_1, &o );

		// non-attached observer
		Observer o2( 3000 );

		// attempt to attach already attached
		TS_ASSERT( link == s.attach( &Observer::on_update_1, &o ) );

		// attempt to detach non-existing observer, different object
		TS_ASSERT( !s.detach( &Observer::on_update_1, &o2 ) );

		// attempt to detach non-existing observer, same object, different function
		TS_ASSERT( !s.detach( &Observer::on_update_2, &o ) );

		// cleanup, remove all observers
		s.clear();
	}


	/// @brief single object multiple connection
	void test_SignalHub_single_to_multi() {
		using namespace sig;

		Subject s;

		Observer o( 5000 );
		s.attach( &Observer::on_update_1, &o );
		s.attach( &Observer::on_update_2, &o );

		Event e;
		e.type = 22;
		e.data = 333;

		s.notify( e );

		TS_ASSERT_EQUALS( 22, o.event1().type );
		TS_ASSERT_EQUALS( 333, o.event1().data );

		TS_ASSERT_EQUALS( 22, o.event2().type );
		TS_ASSERT_EQUALS( 333, o.event2().data );

		// cleanup, remove all observers
		s.clear();
	}


	/// @brief test Link functionality
	void test_SignalHub_link() {
		using namespace sig;
		using utility::signals::Link;

		Link link;
		TS_ASSERT( link.empty() );

		Subject s;

		Observer o( 6000 );
		link = s.attach( &Observer::on_update_1, &o );

		Event e;
		e.type = 55;
		e.data = 789;

		// send first signal
		s.notify( e );

		// invalidate the link
		link.invalidate();
		TS_ASSERT( !link.valid() );

		// attempt to send a second signal
		e.type = 99;
		e.data = 123;
		s.notify( e );

		// signal should not have gone through,
		// obs should be removed
		TS_ASSERT_EQUALS( 55, o.event1().type );
		TS_ASSERT_EQUALS( 789, o.event1().data );
		TS_ASSERT_EQUALS( 0, s.n_obs() );

		// now check invalidation from subject side
		{ // begin scope
			Subject s2;
			link = s2.attach( &Observer::on_update_1, &o );
			TS_ASSERT( link.valid() );
		} // end scope

		// link should be invalid after Subject destroyed
		TS_ASSERT( !link.valid() );
	}

}; // class SignalHubTests

