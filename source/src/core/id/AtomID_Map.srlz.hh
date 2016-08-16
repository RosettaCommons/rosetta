// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/AtomID_Map.srlz.hh
/// @brief  Serialization routines for AtomID_Map
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_id_AtomID_Map_SRLZ_HH
#define INCLUDED_core_id_AtomID_Map_SRLZ_HH

// Unit headers
#include <core/id/AtomID_Map.hh>

// Package headers
#include <core/id/AtomID.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

#ifdef SERIALIZATION
#include <cereal/cereal.hpp>
#endif

namespace core {
namespace id {



#ifdef    SERIALIZATION
/// @brief Serialization routine.  In order for this to successfully compile,
/// the code that is trying to serialize an instance of this class will need
/// to #include <utility/serialization/serialization.hh>, which should already
/// happen, and also #include <utility/vector1.srlz.hh>.
template < class Archive, class T >
void save( Archive & arc, AtomID_Map< T > const & map ) {
	arc( map.default_value() );
	arc( map.size() );
	for ( Size ii = 1; ii <= map.size(); ++ii ) {
		arc( map[ ii ] );
	}
}

/// @brief Deserialization routine; see comments for save, above.
template < class Archive, class T >
void load( Archive & arc, AtomID_Map< T > & map ) {
	T default_val; arc( default_val );
	map.default_value( default_val );
	Size size; arc( size );
	map.resize( size );
	for ( Size ii = 1; ii <= size; ++ii ) {
		arc( map[ ii ] );
	}
}

#endif // SERIALIZATION


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_AtomID_Map_SRLZ_HH
