// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/SingletonBase.hh
/// @brief  A base class for all singltons using CRTP and managing the complexity of safely
///         initializing singltons in a thread-safe way.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_SingletonBase_HH
#define INCLUDED_utility_SingletonBase_HH

#if defined MULTI_THREADED && defined CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif

namespace utility {

/// @brief SingletonBase is meant to serve as a base class for singleton classes in Rosetta
/// handling the initialization of the singleton in a thread-safe way.
///
/// @details The derived class must a) implement a private, static function:
/// T * create_singleton_instance()
/// so that the SingletonBase class can invoke this function, and b) declare the
/// SingletonBase class to be a friend, so that it can invoke this function
/// The .cc file in which the derived singleton must be put will need to include
/// the definitions for the two static data members, instance_ and singleton_mutex_.
template < class T >
class SingletonBase
{
public:
	/// @brief public constructor (the derived class must have a private constructor, of course).
	SingletonBase() {}

	/// @brief Safely instantiate a singleton class in a (possibly)
	/// multithreaded context.
	///
	/// @details In the non-multithreaded case, this simply checks the
	/// singleton's instance member; in the multithreaded case, it
	/// checks the instance member, then it obtains the singleton's
	/// instance-creation mutex, then it checks the instance member
	/// again, to ensure that no other thread has already created the
	/// instance, it creates the instance, and then it releases the mutex.
	static
	T *
	get_instance() {

#ifdef MULTI_THREADED
#ifdef CXX11
		// threadsafe double-checked locking version that uses c++11 interface
		// taken from http://preshing.com/20130930/double-checked-locking-is-fixed-in-cpp11/
		T * local_instance = instance_.load( std::memory_order_relaxed );
		std::atomic_thread_fence( std::memory_order_acquire );
		if ( ! local_instance ) {
			std::lock_guard< std::mutex > lock( singleton_mutex_ );
			local_instance = instance_.load( std::memory_order_relaxed );
			if ( ! local_instance ) {
				local_instance = T::create_singleton_instance();
				instance_.store( local_instance, std::memory_order_relaxed );
				std::atomic_thread_fence( std::memory_order_release );
			}
		}
#else
	  // ok, multithreaded w/o cxx11. Foldit version?
		// Not actually threadsafe!
		// http://www.cs.umd.edu/~pugh/java/memoryModel/DoubleCheckedLocking.html
		if ( ! instance_ ) {
			instance_ = T::create_singleton_instance();
		}
#endif
#else
  	// not multithreaded; standard singleton instantiation logic
		if ( ! instance_ ) {
			instance_ = T::create_singleton_instance();
		}
#endif
		return instance_;
	}

private:
	/// @brief Private, unimplemented copy constructor -- uncopyable.
	SingletonBase( SingletonBase< T > const & );

	/// @brief Private, unimplemented assignment operator -- uncopyable.
	SingletonBase< T > const & operator = ( SingletonBase< T > const & rhs );

private:

#if defined MULTI_THREADED && defined CXX11
	static std::mutex singleton_mutex_;
	static std::atomic< T * > instance_;
#else
	static T * instance_;
#endif

};


} // utility

#endif
