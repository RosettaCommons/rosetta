// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/datacache/ObserverCache.cc
/// @brief  A DataCache storing objects derived from
///         core::pose::datacache::CacheableData.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/pose/datacache/ObserverCache.hh>

// project headers

#include <core/pose/datacache/CacheableObserver.hh>
#include <utility/vector1.hh>
#include <iostream>


namespace core {
namespace pose {
namespace datacache {


/// @brief constructor
/// @param[in] n_types The number of slots for this ObserverCache.
/// @param[in] pose The Pose that will be watched by the Observers in this cache.
ObserverCache::ObserverCache(
	Size const n_slots,
	Pose & pose
) :
	Super( n_slots ),
	pose_( &pose )
{}


/// @brief default destructor
ObserverCache::~ObserverCache() {
	detach();
}


/// @brief copy assignment
ObserverCache & ObserverCache::operator =( ObserverCache const & rval ) {
	if ( this != &rval ) {
		// detach all existing observers
		detach();

		// copy assign (all observers cloned)
		Super::operator =( rval );

		// attach the observers that were already attached in rval
		for ( Size i = 1, ie = rval.data().size(); i <= ie; ++i ) {
			if ( rval.is_attached( i ) ) {
				attach( i );
			}
		}
	}
	return *this;
}


/// @brief clear all the observers
void ObserverCache::clear() {
	detach(); // safety
	Super::clear();
}


/// @brief clear the observer in a selected slot
void ObserverCache::clear( Size const slot ) {
	detach( slot ); // safety
	Super::clear( slot );
}


/// @brief store a copy of the observer in the given slot and attach it to
///  the Pose
/// @param[in] The slot to use.
/// @param[in] observer The Observer to clone() and store.
/// @remarks this function exists to ensure the base class version is
///  overridden
void ObserverCache::set(
	Size const slot,
	CacheableObserverOP observer
)
{
	set( slot, observer, true );
}


/// @brief store a copy of the observer in the given slot
/// @param[in] The slot to use.
/// @param[in] observer The Observer to clone() and store.
/// @param[in] auto_attach Attach the observer to the Pose?
void ObserverCache::set(
	Size const slot,
	CacheableObserverOP observer,
	bool const auto_attach
)
{
	detach( slot ); // safety

	if ( observer.get() != 0 ) {
		data()[ slot ] = observer->clone();

		if ( auto_attach ) {
			attach( slot );
		}
	} else {
		data()[ slot ] = 0;
	}
}


/// @brief is the observer in the slot attached to the Pose?
/// @return true if attached, false if not attached or no observer
///  exists in the slot
bool ObserverCache::is_attached( Size const slot ) const {
	if ( has( slot ) ) {
		return data()[ slot ]->is_attached();
	}

	return false;
}


/// @brief attach all stored observers to the Pose
void ObserverCache::attach() {
	for ( Size i = 1, ie = data().size(); i <= ie; ++i ) {
		attach( i );
	}
}


/// @brief detach all observers from the Pose
void ObserverCache::detach() {
	for ( Size i = 1, ie = data().size(); i <= ie; ++i ) {
		detach( i );
	}
}


/// @brief attach an observer in a particular slot to the Pose
/// @param[in] slot Attach the observer in this slot.
void ObserverCache::attach( Size const slot ) {
	if ( has( slot ) ) {
		data()[ slot ]->attach_to( *pose_ );
	}
}


/// @brief detach an observer in a particular slot to the Pose
/// @param[in] slot Detach the observer in this slot.
void ObserverCache::detach( Size const slot ) {
	if ( has( slot ) ) {
		data()[ slot ]->detach_from();
	}
}


} // namespace datacache
} // namespace pose
} // namespace core
