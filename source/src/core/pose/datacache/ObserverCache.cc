// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace datacache {


/// @brief constructor
/// @param[in] n_types The number of slots for this ObserverCache.
/// @param[in] pose The Pose that will be watched by the Observers in this cache.
ObserverCache::ObserverCache(
	Size n_slots,
	Pose & pose
) :
	Super( n_slots ),
	attached_( n_slots, false ),
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
void ObserverCache::clear( Size slot ) {
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
	Size slot,
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
	Size slot,
	CacheableObserverOP observer,
	bool auto_attach
)
{
	detach( slot ); // safety

	if ( observer.get() != nullptr ) {
		data()[ slot ] = observer->clone();

		if ( auto_attach ) {
			attach( slot );
		}
	} else {
		data()[ slot ] = nullptr;
	}
}


/// @brief is the observer in the slot attached to the Pose?
/// @return true if attached, false if not attached or no observer
///  exists in the slot
bool ObserverCache::is_attached( Size slot ) const {
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
void ObserverCache::attach( Size slot ) {
	if ( has( slot ) ) {
		data()[ slot ]->attach_to( *pose_ );
		attached_[ slot ] = true;
	}
}


/// @brief detach an observer in a particular slot to the Pose
/// @param[in] slot Detach the observer in this slot.
void ObserverCache::detach( Size slot ) {
	if ( has( slot ) ) {
		data()[ slot ]->detach_from();
		attached_[ slot ] = false;
	}
}


} // namespace datacache
} // namespace pose
} // namespace core

#ifdef    SERIALIZATION

/// @details The default constructor is provided only for deserialization
core::pose::datacache::ObserverCache::ObserverCache() :
	Super(),
	pose_( 0 )
{}

/// @details Following deserialization, attach a Pose to this observer (or the observer
/// to the pose, really) and then make sure that all of the CacheableObservers that
/// were attached prior to serialization become reattached.
void
core::pose::datacache::ObserverCache::attach_pose( Pose & pose )
{
	debug_assert( pose_ == 0 );
	pose_ = &pose;
	for ( Size ii = 1; ii <= attached_.size(); ++ii ) {
		if ( attached_[ ii ] ) {
			attach( ii );
		}
	}
}


/// @brief Serialization method
/// @details This class does not serialize the raw pointer to the pose that it is
/// observing; it will be the responsibility of the Pose that owns this cache to
/// initialize the deserialized %ObserverCache with a pointer to itself when that
/// pose itself is being deserialized.
template< class Archive >
void
core::pose::datacache::ObserverCache::save( Archive & arc ) const {
	arc( data() );
	arc( attached_ );
	// The raw pointer to the pose_ is not savable
	// EXEMPT pose_
}

/// @brief Deserialization method; it will not attempt to deserialize the raw pointer
/// to the Pose that it's meant to observer; rather, responsibility for initializing
/// that pointer falls on the Pose that holds this %ObserverCache
template< class Archive >
void
core::pose::datacache::ObserverCache::load( Archive & arc ) {
	arc( data() );
	arc( attached_ );
	// The raw pointer to the pose_ is not savable
	// EXEMPT pose_
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::datacache::ObserverCache );
CEREAL_REGISTER_TYPE( core::pose::datacache::ObserverCache )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_datacache_ObserverCache )
#endif // SERIALIZATION
