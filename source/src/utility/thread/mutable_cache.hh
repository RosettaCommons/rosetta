// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/thread/mutable_cache.hh
/// @brief  A key:value cache intended to be used for (threadsafe) mutable caching of values.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_utility_thread_mutable_cache_HH
#define INCLUDED_utility_thread_mutable_cache_HH

// Unit headers

// Utility headers

// Boost headers
#include <functional>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#include <thread>
#endif

// C++ headers
#include <map>

namespace utility {
namespace thread {

template< class K, class V >
class MutableCache {
public:

	MutableCache() = default;
	~MutableCache() = default;
	MutableCache( MutableCache const & src );
	MutableCache( MutableCache && src );
	MutableCache & operator=( MutableCache const & src );
	MutableCache & operator=( MutableCache && src );

public:

	bool has( K const & key ) const;

	V const & get( K const & key ) const;

	/// @brief Add the key:value pair to the map, but only if the key isn't already present.
	void add_if_missing( K const & key, V const & value );

	/// @brief Add a key:value pair to the map if the key isn't present.
	/// The creator function will be called to generate the value
	void add_if_missing( K const & key, std::function< V() > creator );

	/// @brief Add a key:value pair to the map if the key isn't present.
	/// The creator function will be called to generate the value
	void add_if_missing( K const & key, std::function< V(K const &) > creator );

private:
#ifdef MULTI_THREADED
	mutable ReadWriteMutex mutex_;
#endif

	std::map< K, V > data_;
};

// Out-of-line definitions
template< class K, class V >
MutableCache<K,V>::MutableCache( MutableCache const & src ) {
#ifdef MULTI_THREADED
	ReadLockGuard rlg( src.mutex_ );
#endif
	data_ = src.data_;
}

template< class K, class V >
MutableCache<K,V>::MutableCache( MutableCache && src ):
	data_( std::move( src.data_ ) )
{}

template< class K, class V >
MutableCache<K,V> &
MutableCache<K,V>::operator=( MutableCache const & src ) {
	if ( this == &src ) return *this; // Check self assignment to avoid deadlock
#ifdef MULTI_THREADED
	PairedReadLockWriteLockGuard rwlg( src.mutex_, this->mutex_ );
#endif
	data_ = src.data_;
	return *this;
}

template< class K, class V >
MutableCache<K,V> &
MutableCache<K,V>::operator=( MutableCache && src ) {
#ifdef MULTI_THREADED
	WriteLockGuard wlg( this->mutex_ );
	// rvalue ref implies no one else has access to it.
#endif
	data_ = std::move( src.data_ );
	return *this;
}

template< class K, class V >
bool
MutableCache<K,V>::has( K const & key ) const {
#ifdef MULTI_THREADED
	ReadLockGuard rlg( mutex_ );
#endif
	return data_.count( key ) > 0;
}

template< class K, class V >
V const &
MutableCache<K,V>::get( K const & key ) const {
#ifdef MULTI_THREADED
	ReadLockGuard rlg( mutex_ );
#endif
	return data_.at(key);
}

template< class K, class V >
void
MutableCache<K,V>::add_if_missing( K const & key, V const & value ) {
#ifdef MULTI_THREADED
	WriteLockGuard rlg( mutex_ );
#endif
	if ( data_.count( key ) == 0 ) {
		data_[key] = value;
	}
}

template< class K, class V >
void
MutableCache<K,V>::add_if_missing( K const & key, std::function< V() > creator ) {
#ifdef MULTI_THREADED
	WriteLockGuard rlg( mutex_ );
#endif
	if ( data_.count( key ) == 0 ) {
		data_[key] = creator();
	}
}

template< class K, class V >
void
MutableCache<K,V>::add_if_missing( K const & key, std::function< V(K const &) > creator ) {
#ifdef MULTI_THREADED
	WriteLockGuard rlg( mutex_ );
#endif
	if ( data_.count( key ) == 0 ) {
		data_[key] = creator(key);
	}
}


} //namespace thread
} //namespace utility

#endif
