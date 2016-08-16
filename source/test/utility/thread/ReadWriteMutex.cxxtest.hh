// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/thread/ReadWriteMutex.cxxtest.hh
/// @brief  test suite for the ReadWriteMutex class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Package headers
#include <utility/LexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.tmpl.hh>
#include <cxxtest/TestSuite.h>

/// C++ headers
#include <iostream>

#include <utility/thread/ReadWriteMutex.hh>
#include <utility/string_util.hh>

#ifdef MULTI_THREADED
#ifdef CXX11

#include <thread>
#include <chrono>


void
test_two_threads_read_together_thread_1(
	utility::thread::ReadWriteMutex & wrm,
	utility::vector1< std::mutex > & mutexes,
	utility::vector1< std::string > & read_sequence
)
{
	mutexes[1].lock();
	{ // scope
	utility::thread::ReadLockGuard rlg( wrm );
	read_sequence.push_back( "Thread 1 read" );
	mutexes[2].unlock();
  mutexes[3].lock();
  read_sequence.push_back( "Thread 1 about to release ReadLockGuard" );
	}
  mutexes[4].unlock();
}

void
test_two_threads_read_together_thread_2(
	utility::thread::ReadWriteMutex & wrm,
	utility::vector1< std::mutex > & mutexes,
	utility::vector1< std::string > & read_sequence
)
{
	mutexes[2].lock();
	{ // scope
	utility::thread::ReadLockGuard rlg( wrm );
	read_sequence.push_back( "Thread 2 read" );
	mutexes[3].unlock();
  mutexes[4].lock();
  read_sequence.push_back( "Thread 2 about to release ReadLockGuard" );
	}
}

void
test_thread_waits_for_reading_to_finish_to_write_t1(
	utility::thread::ReadWriteMutex & wrm,
	utility::vector1< std::mutex > & mutexes,
	utility::vector1< std::string > & read_sequence
)
{
	mutexes[1].lock();
	read_sequence.push_back( "Thread 1: initial read counter = " + utility::to_string( wrm.read_counter() ) );
	{
	utility::thread::ReadLockGuard rlg( wrm );
	read_sequence.push_back( "Thread 1 read begin (read counter now = " + utility::to_string( wrm.read_counter() ) + ")" );
	mutexes[2].unlock();
	mutexes[3].lock();
	std::chrono::milliseconds dura( 50 );
	std::this_thread::sleep_for( dura );
	read_sequence.push_back( "Thread 1 about to release ReadLockGuard" );
	}
}

void
test_thread_waits_for_reading_to_finish_to_write_t2(
	utility::thread::ReadWriteMutex & wrm,
	utility::vector1< std::mutex > & mutexes,
	utility::vector1< std::string > & read_sequence
)
{
	mutexes[2].lock();
	{
	read_sequence.push_back( "Thread 2 try and aquire write lock (read counter now = " + utility::to_string( wrm.read_counter() ) + ")" );
	mutexes[3].unlock();
	utility::thread::WriteLockGuard wlg( wrm );
	read_sequence.push_back( "Thread 2 aquired write lock" );
	}
}

void
test_write_lock_prevents_reading_t1(
	utility::thread::ReadWriteMutex & wrm,
	utility::vector1< std::mutex > & mutexes,
	utility::vector1< std::string > & read_sequence
)
{
	mutexes[1].lock();
	{
		read_sequence.push_back( "Thread 1 about to obtain write lock" );
		utility::thread::WriteLockGuard guard( wrm );
		mutexes[2].unlock();
		std::chrono::milliseconds dura( 50 );
		std::this_thread::sleep_for( dura );
		read_sequence.push_back( "Thread 1 waking up, about to release write lock" );
	}
}

void
test_write_lock_prevents_reading_t2(
	utility::thread::ReadWriteMutex & wrm,
	utility::vector1< std::mutex > & mutexes,
	utility::vector1< std::string > & read_sequence
)
{
	mutexes[2].lock();
	{
		read_sequence.push_back( "Thread 2 about to obtain read lock" );
		utility::thread::ReadLockGuard guard( wrm );
		read_sequence.push_back( "Thread 2 finally obtained read lock" );
	}
}


#endif
#endif

class ReadWriteMutexTests : public CxxTest::TestSuite {
public:
	void test_empty() {
		// appease compiler in case we didn't build with the extras=cxx11thread flag
		TS_ASSERT( true );
	}


public:

	/// @brief ensure that two threads can read from the
	void test_two_threads_read_together() {
#ifdef MULTI_THREADED
#ifdef CXX11

		utility::thread::ReadWriteMutex rwm;

		utility::vector1< std::string > read_sequence;
		utility::vector1< std::mutex > mutexes( 4 );
		for ( int ii = 1; ii <= 4; ++ii ) mutexes[ ii ].lock();

		std::thread t1( test_two_threads_read_together_thread_1, std::ref(rwm), std::ref(mutexes), std::ref(read_sequence) );
		std::thread t2( test_two_threads_read_together_thread_2, std::ref(rwm), std::ref(mutexes), std::ref(read_sequence) );

		mutexes[1].unlock();
		t1.join();
		t2.join();

		TS_ASSERT( read_sequence.size() == 4 );
		TS_ASSERT( read_sequence[1] == "Thread 1 read" );
		TS_ASSERT( read_sequence[2] == "Thread 2 read" );
		TS_ASSERT( read_sequence[3] == "Thread 1 about to release ReadLockGuard" );
		TS_ASSERT( read_sequence[4] == "Thread 2 about to release ReadLockGuard" );
		//for ( unsigned int ii = 1; ii <= read_sequence.size(); ++ii ) {
		//	std::cout << "ii = " << ii << " " << read_sequence[ii] << std::endl;
		//}
#endif
#endif

	}

	void test_thread_waits_for_reading_to_finish_to_write() {
#ifdef MULTI_THREADED
#ifdef CXX11

		utility::thread::ReadWriteMutex rwm;

		utility::vector1< std::string > read_sequence;
		utility::vector1< std::mutex > mutexes( 4 );
		for ( int ii = 1; ii <= 4; ++ii ) mutexes[ ii ].lock();

		std::thread t1( test_thread_waits_for_reading_to_finish_to_write_t1, std::ref(rwm), std::ref(mutexes), std::ref(read_sequence) );
		std::thread t2( test_thread_waits_for_reading_to_finish_to_write_t2, std::ref(rwm), std::ref(mutexes), std::ref(read_sequence) );
		mutexes[1].unlock();
		t1.join();
		t2.join();

		TS_ASSERT( read_sequence.size() == 5 );
		TS_ASSERT( read_sequence[1] == "Thread 1: initial read counter = 0" );
		TS_ASSERT( read_sequence[2] == "Thread 1 read begin (read counter now = 1)" );
		TS_ASSERT( read_sequence[3] == "Thread 2 try and aquire write lock (read counter now = 1)" );
		TS_ASSERT( read_sequence[4] == "Thread 1 about to release ReadLockGuard" );
		TS_ASSERT( read_sequence[5] == "Thread 2 aquired write lock" );
		//for ( unsigned int ii = 1; ii <= read_sequence.size(); ++ii ) {
		//	std::cout << "ii = " << ii << " " << read_sequence[ii] << std::endl;
		//}
#endif
#endif

	}

	void test_write_lock_prevents_reading() {
#ifdef MULTI_THREADED
#ifdef CXX11

		utility::thread::ReadWriteMutex rwm;

		utility::vector1< std::string > read_sequence;
		utility::vector1< std::mutex > mutexes( 2 );
		for ( int ii = 1; ii <= 2; ++ii ) mutexes[ ii ].lock();

		std::thread t1( test_write_lock_prevents_reading_t1, std::ref(rwm), std::ref(mutexes), std::ref(read_sequence) );
		std::thread t2( test_write_lock_prevents_reading_t2, std::ref(rwm), std::ref(mutexes), std::ref(read_sequence) );
		mutexes[1].unlock();
		t1.join();
		t2.join();

		TS_ASSERT( read_sequence.size() == 4 );
		TS_ASSERT( read_sequence[1] == "Thread 1 about to obtain write lock" );
		TS_ASSERT( read_sequence[2] == "Thread 2 about to obtain read lock" );
		TS_ASSERT( read_sequence[3] == "Thread 1 waking up, about to release write lock" );
		TS_ASSERT( read_sequence[4] == "Thread 2 finally obtained read lock" );
		//for ( unsigned int ii = 1; ii <= read_sequence.size(); ++ii ) {
		//	std::cout << "ii = " << ii << " " << read_sequence[ii] << std::endl;
		//}
#endif
#endif

	}

};
