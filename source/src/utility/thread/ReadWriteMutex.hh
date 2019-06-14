// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/thread/ReadWriteMutex.hh
/// @brief  Classes to manage data that can be read by multiple threads and written to by only one thread
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Modified by Vikram K. Mulligan (vmulligan@flatironinstitute.org) to add additional options for
/// conditionally locking mutexes or for locking multiple mutexes simultaneously without deadlock.

#ifndef INCLUDED_utility_thread_ReadWriteMutex_HH
#define INCLUDED_utility_thread_ReadWriteMutex_HH

#ifdef MULTI_THREADED

// Unit headers
#include <utility/thread/ReadWriteMutex.fwd.hh>

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
	/// Attempting to obtain two read locks in the same thread
	/// (e.g. in one function and then again in a called function) can result in deadlock.
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
	std::mutex read_lock_;
	std::atomic< int > read_counter_;
	std::mutex write_lock_;
	std::condition_variable_any not_being_read_;
};

/// @brief A class for RAII-based read-locking of ReadWriteMutexes.
/// @author Andrew Leaver-Fay.
/// @author Modified by Vikram K. Mulligan (vmulligan@flatironinstitute.org) to add do_nothing option.
class ReadLockGuard
{
public:
	/// @brief Block until you obtain read permission
	/// @details Do not instantiate if the current thread has any other read or write locks for this mutex.
	/// @note If do_nothing is true, this does nothing.  This is useful in cases in which you might want to specify
	/// that the mutex has already been locked.
	ReadLockGuard( ReadWriteMutex & rwm, bool const do_nothing = false );

	/// @brief Default constructor is explicitly deleted.
	ReadLockGuard() = delete;

	/// @brief Copy constructor is explicitly deleted.
	ReadLockGuard( ReadLockGuard const & ) = delete;

	/// @brief Assignment operator is explicitly deleted.
	ReadLockGuard & operator= ( ReadLockGuard const & ) = delete;

	/// @brief Invoke release_read_lock() on the RealWriteLock used to construct this instance
	~ReadLockGuard();
private:
	ReadWriteMutex & rwm_;
	bool did_nothing_ = false;
};

/// @brief A class for RAII-based write-locking of ReadWriteMutexes.
/// @author Andrew Leaver-Fay.
/// @author Modified by Vikram K. Mulligan (vmulligan@flatironinstitute.org) to add do_nothing option.
class WriteLockGuard
{
public:
	/// @brief Block until you obtain write permission
	/// @details Do not instantiate if the current thread has any other read or write locks for this mutex.
	/// @note If do_nothing is true, this does nothing.  This is useful in cases in which you might want to specify
	/// that the mutex has already been locked.
	WriteLockGuard( ReadWriteMutex & rwm, bool const do_nothing = false );

	/// @brief Default constructor is explicitly deleted.
	WriteLockGuard() = delete;

	/// @brief Copy constructor is explicitly deleted.
	WriteLockGuard( WriteLockGuard const & ) = delete;

	/// @brief Assignment operator is explicitly deleted.
	WriteLockGuard & operator= ( WriteLockGuard const & ) = delete;

	/// @brief Invoke release_write_lock() on the ReadWriteMutex used to construct this instance
	~WriteLockGuard();

private:
	ReadWriteMutex & rwm_;
	bool did_nothing_ = false;

};

/// @brief A class for RAII-based read- or write-locking of ReadWriteMutexes.
/// @details The choice of whether to read-lock or write-lock must be made at object instantiation time.#pragma endregion
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class ReadOrWriteLockGuard
{
public:
	/// @brief Constructor.  If read_mode is true, this obtains a read lock on the mutex.  If it's false,
	/// it obtains a write-lock.  If do_nothing is true, the behaviour is overridden and the lock guard
	/// does nothing.
	ReadOrWriteLockGuard( ReadWriteMutex & rwm, bool const read_mode, bool const do_nothing = false );

	/// @brief Default constructor is explicitly deleted.
	ReadOrWriteLockGuard() = delete;

	/// @brief Copy constructor is explicitly deleted.
	ReadOrWriteLockGuard( ReadOrWriteLockGuard const & ) = delete;

	/// @brief Assignment operator is explicitly deleted.
	ReadOrWriteLockGuard & operator= ( ReadOrWriteLockGuard const & ) = delete;

	/// @brief Default destructor.
	~ReadOrWriteLockGuard() = default;

private:
	utility::pointer::shared_ptr< ReadLockGuard > read_guard_ = nullptr;
	utility::pointer::shared_ptr< WriteLockGuard > write_guard_ = nullptr;
};


/// @brief A RAII lock guard class that simultaneously read-locks one mutex and write-locks another, while avoiding deadlock.
/// @details The choice of which mutex to lock first depends on the addresses of the two mutexes.  The lower memory address is
/// always locked first, preventing cases in which thread 1 has locked A and is waiting for a lock on B while thread 2 has locked
/// B and is waiting on a lock for A.  If this class is used, both thread 1 and thread 2 will choose the same mutex to lock first.
/// @note This is useful in cases in which a thread might want to read-lock object A and write-lock object B to copy data from A
/// to B, if there is a risk that another thread might want to copy data from B to A.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class PairedReadLockWriteLockGuard
{
public:
	/// @brief Constructor.  The first mutex gets a read lock, and the second gets a write lock.  If do_nothing is true,
	/// the behaviour is overridden and the lock guard does nothing.
	PairedReadLockWriteLockGuard( ReadWriteMutex & mutex_to_read_lock, ReadWriteMutex & mutex_to_write_lock, bool const do_nothing = false );

	/// @brief Default constructor is explicitly deleted.
	PairedReadLockWriteLockGuard() = delete;

	/// @brief Copy constructor is explicitly deleted.
	PairedReadLockWriteLockGuard( PairedReadLockWriteLockGuard const & ) = delete;

	/// @brief Assignment operator is explicitly deleted.
	PairedReadLockWriteLockGuard & operator= ( PairedReadLockWriteLockGuard const & ) = delete;

	/// @brief Default destructor.
	~PairedReadLockWriteLockGuard() = default;

private:

	utility::pointer::shared_ptr< ReadOrWriteLockGuard > first_guard_ = nullptr;
	utility::pointer::shared_ptr< ReadOrWriteLockGuard > second_guard_ = nullptr;

};

}
}

#endif

#endif
