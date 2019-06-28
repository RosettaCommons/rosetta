// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/thread/ReadWriteMutex.cc
/// @brief  Classes to manage data that can be read by multiple threads and written to by only one thread
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Modified by Vikram K. Mulligan (vmulligan@flatironinstitute.org) to add additional options for
/// conditionally locking mutexes or for locking multiple mutexes simultaneously without deadlock.

#ifdef MULTI_THREADED

// Unit headers
#include <utility/thread/ReadWriteMutex.hh>

// Error messages
#include <utility/exit.hh>

// C++ Headers
#include <cassert>

namespace utility {
namespace thread {

ReadWriteMutex::ReadWriteMutex() :
	read_counter_( 0 )
{}

void ReadWriteMutex::obtain_read_lock()
{
	std::lock_guard< std::mutex > lock( read_lock_ );
	++read_counter_;
}

void ReadWriteMutex::release_read_lock()
{
	assert( read_counter_ > 0 );
	--read_counter_;
	if ( read_counter_ == 0 ) {
		std::lock_guard< std::mutex > lock( write_lock_ );
		not_being_read_.notify_one();
	}
}

void ReadWriteMutex::obtain_write_lock()
{
	std::unique_lock< std::mutex > prevent_more_reads( read_lock_ );
	std::unique_lock< std::mutex > no_more_readers( write_lock_ );

	not_being_read_.wait( no_more_readers, [ this ]() { return this->read_counter_.load() == 0; } );

	// Ok, we now have both the read and the write locks and the read counter is at 0.
	// Release the mutexes without unlocking them and exit this function

	prevent_more_reads.release();
	no_more_readers.release();
}

void ReadWriteMutex::release_write_lock()
{
	read_lock_.unlock();
	write_lock_.unlock();
}


ReadLockGuard::ReadLockGuard( ReadWriteMutex & rwm, bool const do_nothing /*= false*/ ) :
	rwm_( rwm ),
	did_nothing_( do_nothing )
{
	if ( !do_nothing ) {
		rwm_.obtain_read_lock();
	}
}

ReadLockGuard::~ReadLockGuard()
{
	if ( !did_nothing_ ) {
		rwm_.release_read_lock();
	}
}

WriteLockGuard::WriteLockGuard( ReadWriteMutex & rwm, bool const do_nothing /*= false*/ ) :
	rwm_( rwm ),
	did_nothing_( do_nothing )
{
	if ( !do_nothing ) {
		rwm_.obtain_write_lock();
	}
}

WriteLockGuard::~WriteLockGuard()
{
	if ( !did_nothing_ ) {
		rwm_.release_write_lock();
	}
}

/// @brief Constructor.  If read_mode is true, this obtains a read lock on the mutex.  If it's false,
/// it obtains a write-lock.  If do_nothing is true, the behaviour is overridden and the lock guard
/// does nothing.
ReadOrWriteLockGuard::ReadOrWriteLockGuard(
	ReadWriteMutex & rwm,
	bool const read_mode,
	bool const do_nothing /*=false*/
) {
	if ( !do_nothing ) {
		if ( read_mode ) {
			read_guard_ = utility::pointer::make_shared< ReadLockGuard >( rwm );
		} else {
			write_guard_ = utility::pointer::make_shared< WriteLockGuard >( rwm );
		}
	}
}

/// @brief Constructor.  The first mutex gets a read lock, and the second gets a write lock.  If do_nothing is true,
/// the behaviour is overridden and the lock guard does nothing.
PairedReadLockWriteLockGuard::PairedReadLockWriteLockGuard(
	ReadWriteMutex & mutex_to_read_lock,
	ReadWriteMutex & mutex_to_write_lock,
	bool const do_nothing /*= false*/
) {
	if ( !do_nothing ) {
		runtime_assert_string_msg( &mutex_to_read_lock != &mutex_to_write_lock, "Error in utility::thread::PairedReadLockWriteLockGuard(): Two copies of the same mutex were passed to this function, which will inevitibly cause deadlock!" );

		//Objective criterion for chosing which mutex to lock first (ensuring that all threads try to lock the same mutex first):
		bool const read_mutex_is_greater( &mutex_to_read_lock > &mutex_to_write_lock );
		//The following may look odd, but it's necessary to ensure that deadlock doesn't occur if ever
		//two different threads do a=b and b=a simultaneously.  Basically, we're ensuring that both threads
		//try to lock the same mutex first, so that there isn't a deadlock scenario in which X hold mutex A
		//and is waiting to acquire mutex B while Y holds mutex B and is waiting to acquire mutex A.  This way,
		//either X or Y will get mutex A, then will proceed to acquire mutex B and the other will block until
		//the first is done.
		first_guard_ = utility::pointer::make_shared< ReadOrWriteLockGuard >( read_mutex_is_greater ? mutex_to_write_lock : mutex_to_read_lock, !read_mutex_is_greater );
		second_guard_ = utility::pointer::make_shared< ReadOrWriteLockGuard >( read_mutex_is_greater ? mutex_to_read_lock : mutex_to_write_lock, read_mutex_is_greater );
	}
}

}
}

#endif

