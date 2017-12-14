// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/datacache/CacheableObserver.cc
/// @brief  Base class for Pose/Conformation observers that are stored in
///         a Pose's DataCache.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/pose/datacache/CacheableObserver.hh>
#include <core/pose/datacache/ObserverCache.hh>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace datacache {


/// @brief default constructor
CacheableObserver::CacheableObserver() :
	Super()
{}


/// @brief copy constructor
CacheableObserver::CacheableObserver( CacheableObserver const & /*rval*/ ) = default;


/// @brief default destructor
/// @warning Derived classes must remember to detach on destruction!
CacheableObserver::~CacheableObserver() = default;


/// @brief copy assignment
CacheableObserver &
CacheableObserver::operator =( CacheableObserver const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
	}

	return *this;
}



/// @brief attach to Pose/Conformation
///  Derived classes do not overload this method -- see attach_impl()
///  instead.
void CacheableObserver::attach_to( Pose & pose ) {
	detach_from();
	attach_impl( pose );
}


/// @brief detach from Pose/Conformation
/// @remarks Derived classes do not overload this method -- see
///  detach_impl() instead.
void CacheableObserver::detach_from() {
	detach_impl();
}


} // namespace datacache
} // namespace pose
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::datacache::CacheableObserver::save( Archive & ) const {
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::datacache::CacheableObserver::load( Archive & ) {
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::datacache::CacheableObserver );
CEREAL_REGISTER_TYPE( core::pose::datacache::CacheableObserver )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_datacache_CacheableObserver )
#endif // SERIALIZATION
