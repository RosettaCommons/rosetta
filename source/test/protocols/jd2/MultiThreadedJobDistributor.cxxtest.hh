// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MultiThreadedJobDistributor.cxxtest.hh
/// @brief  test suite for protocols::jd2::MultiThreadedJobDistributor
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
// #include <protocols/jd2/MultiThreadedJobDistributor.hh>

// Utility headers
#include <utility/thread/backwards_thread_local.hh>

// C++ headers
#include <string>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 Headers
#include <mutex>
#include <thread>

#endif
#endif

// using namespace protocols::jd2;

class Dummy {
public:
	Dummy() { my_index_ = ++n_dummies; }
	int index() const { return my_index_; }
private:
	int my_index_;
	static int n_dummies;
};

int Dummy::n_dummies( 0 );

Dummy *
get_dummy() {
	static THREAD_LOCAL Dummy * dummy( 0 );
	if ( dummy == 0 ) {
		dummy = new Dummy;
	}
	return dummy;
}

class ThreadRunner {
public:
	ThreadRunner( int target ) {
		target_ = target;
	}

	void go() {
		Dummy * my_dummy = get_dummy();
		//std::cout << "Thread runner: " << my_dummy->index() << std::endl;
		TS_ASSERT( my_dummy->index() == target_ );
	}

private:
	int target_;

};

class MultiThreadedJobDistributorTests : public CxxTest::TestSuite {
public:

	void setUp() {
		protocols_init();
	}

	// This test is just for the sake of understanding how thread_local works
	void test_dummy_in_tss_one_thread() {
		Dummy * the_dummy = get_dummy();
		TS_ASSERT_EQUALS( the_dummy->index(), 1 );
	}

	// This test is just for the sake of understanding how thread_local works
	void test_tss_in_threads() {
		TS_ASSERT( true );
#ifdef CXX11
#ifdef MULTITHREADED
		ThreadRunner runner( 2 );
		std::thread mythread( &ThreadRunner::go, runner );
		mythread.join();
#endif
#endif
	}

};
