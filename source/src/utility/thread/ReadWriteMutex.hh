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

/// @file   utility/thread/ReadWriteMutex.hh
/// @brief  Classes to manage data that can be read by multiple threads and written to by only one thread
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_thread_ReadWriteMutex_HH
#define INCLUDED_utility_thread_ReadWriteMutex_HH

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 Headers
#include <atomic>
#include <condition_variable>
#include <mutex>

namespace utility {
namespace thread {

class ReadWriteMutex
{
public:
	ReadWriteMutex();

	/// @brief Block until permission to read from a given data structure is granted.
	/// Until this thread has invoked "release_read_lock," there is a guarantee that
	/// no other thread will be able to obtain a write lock.  Multiple threads may
	/// obtain a read lock simultaneously.
	void obtain_read_lock();

	/// @brief After having obtained a read lock, release it.  Never release a read lock
	/// if you have not obtained it first.
	void release_read_lock();

	/// @brief Block until permission to write to a given data structure is granted.
	/// This requires that all threads that have obtained read locks have released them.
	/// No other thread will be able to obtain a read lock or a write lock until this
	/// thread has released its write lock.  Note: do not try to obtain a write lock
	/// if a thread had previously obtained a read lock without first releasing it,
	/// or you will hit a deadlock.
	void obtain_write_lock();

	/// @brief After having obtained the write lock, release it.  Never release a read lock
	/// if you have not obtained it first.
	void release_write_lock();

	int read_counter() const { return read_counter_.load(); }

private:
	std::recursive_mutex read_lock_;
	std::atomic< int > read_counter_;
	std::recursive_mutex write_lock_;
	std::condition_variable_any not_being_read_;
};

class ReadLockGuard
{
public:
	/// @brief Block until you obtain read permission
	ReadLockGuard( ReadWriteMutex & rwm );

	/// @brief Invoke release_read_lock() on the RealWriteLock used to construct this instance
	~ReadLockGuard();
private:
	ReadWriteMutex & rwm_;
};

class WriteLockGuard
{
public:
	/// @brief Block until you obtain write permission
	WriteLockGuard( ReadWriteMutex & rwm );

	/// @brief Invoke release_write_lock() on the ReadWriteMutex used to construct this instance
	~WriteLockGuard();

private:
	ReadWriteMutex & rwm_;

};

}
}

#endif
#endif

#endif
