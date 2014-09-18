// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   utility/thread/ReadWriteMutex.cc
/// @brief  Classes to manage data that can be read by multiple threads and written to by only one thread
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef MULTI_THREADED
#ifdef CXX11

// Unit headers
#include <utility/thread/ReadWriteMutex.hh>


// C++ Headers
#include <cassert>

namespace utility {
namespace thread {

ReadWriteMutex::ReadWriteMutex() :
	read_counter_( 0 )
{}

void ReadWriteMutex::obtain_read_lock()
{
	std::lock_guard< std::recursive_mutex > lock( read_lock_ );
	++read_counter_;
}

void ReadWriteMutex::release_read_lock()
{
	assert( read_counter_ > 0 );
	--read_counter_;
	if ( read_counter_ == 0 ) {
		not_being_read_.notify_one();
	}
}

void ReadWriteMutex::obtain_write_lock()
{
	std::unique_lock< std::recursive_mutex > prevent_more_reads( read_lock_ );
	std::unique_lock< std::recursive_mutex > no_more_readers( write_lock_ );

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


ReadLockGuard::ReadLockGuard( ReadWriteMutex & rwm ) :
	rwm_( rwm )
{
	rwm_.obtain_read_lock();
}

ReadLockGuard::~ReadLockGuard()
{
	rwm_.release_read_lock();
}

WriteLockGuard::WriteLockGuard( ReadWriteMutex & rwm ) :
	rwm_( rwm )
{
	rwm_.obtain_write_lock();
}

WriteLockGuard::~WriteLockGuard()
{
	rwm_.release_write_lock();
}

}
}

#endif
#endif

