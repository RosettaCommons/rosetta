// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/// @file ResidueGraphTypes.srlz.hh
///
/// @brief Helper functions for serializing objects in ResidueGraphTypes.hh and associated
///
/// @author Rocco Moretti (rmorettiase@gmail.com)
////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_ResidueGraphTypes_srlz_hh
#define INCLUDED_core_chemical_ResidueGraphTypes_srlz_hh

#ifdef    SERIALIZATION

// Unit headers
#include <core/chemical/ResidueGraphTypes.hh>

#include <utility/vector1.srlz.hh>

#include <cstdint> // for uintptr_t
#include <cereal/cereal.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#include <cereal/types/map.hpp>

namespace core {
namespace chemical {

// uintptr_t is a C++11 type which should be large enough to hold any pointer (what the VD is in practice)
// It's listed as "optional" though, which may mean that you need a platform/compiler-specific typedef
template< class Archive >
void SERIALIZE_VD( Archive & arc, core::chemical::VD var, std::string const & name ) {
	uintptr_t converted( reinterpret_cast<uintptr_t>(var) );
	if ( var == boost::graph_traits<ResidueGraph>::null_vertex() ) {
		converted = 0; // Want to make sure that null_verticies are correctly translated on the odd case the null_vertex isn't consistent.
	}
	//std::cout << "OUT: " << converted << std::endl;
	arc( ::cereal::make_nvp( name, converted ) );
}
template< class Archive >
void SERIALIZE_VD( Archive & arc, core::chemical::VD var ) {
	uintptr_t converted( reinterpret_cast<uintptr_t>(var) );
	if ( var == boost::graph_traits<ResidueGraph>::null_vertex() ) {
		converted = 0; // Want to make sure that null_verticies are correctly translated.
	}
	//std::cout << "OUT: " << converted << std::endl;
	arc( converted );
}
template< class Archive >
void DESERIALIZE_VD( Archive & arc, core::chemical::VD & var ) {
	uintptr_t converted; arc( converted );
	//std::cout << "IN: " << converted << std::endl;
	if ( converted == 0 ) {
		var = boost::graph_traits<ResidueGraph>::null_vertex();
	} else {
		var = reinterpret_cast< core::chemical::VD >(converted);
	}
}

template< class Archive >
void DESERIALIZE_VD( Archive & arc, core::chemical::VD & var, std::map< VD, VD > const & translation  ) {
	uintptr_t converted; arc( converted );
	//std::cout << "IN: " << converted << std::endl;
	if ( converted == 0 ) {
		var = boost::graph_traits<ResidueGraph>::null_vertex(); // Don't bother to convert the null_vertex.
	} else {
		VD back_convert( reinterpret_cast< core::chemical::VD >(converted) );
		debug_assert( translation.count( back_convert ) == 1 );
		var = translation.at( back_convert );
	}
}

template< class Archive >
void SERIALIZE_VD_VD_MAP( Archive & arc, std::map< VD, VD > const & map ) {
	arc( map.size() );
	for ( std::pair< VD, VD > item : map ) {
		SERIALIZE_VD( arc, item.first );
		SERIALIZE_VD( arc, item.second );
	}
}

template< class Archive >
void DESERIALIZE_VD_VD_MAP( Archive & arc, std::map< VD, VD > & map, std::map< VD, VD > const & translation ) {
	map.clear();
	core::Size nentries; arc( nentries );
	for ( core::Size ii(1); ii <= nentries; ++ii ) {
		VD key, value;
		DESERIALIZE_VD( arc, key, translation );
		DESERIALIZE_VD( arc, value, translation );
		map[ key ] = value;
	}
}

template< class Archive >
void SERIALIZE_VD_VECTOR( Archive & arc, utility::vector1< VD > const & vector ) {
	arc( vector.size() );
	for ( VD item : vector ) {
		SERIALIZE_VD( arc, item );
	}
}

template< class Archive >
void DESERIALIZE_VD_VECTOR( Archive & arc, utility::vector1< VD > & vector, std::map< VD, VD > const & translation ) {
	vector.clear();
	core::Size nentries; arc( nentries );
	for ( core::Size ii(1); ii <= nentries; ++ii ) {
		VD value;
		DESERIALIZE_VD( arc, value, translation );
		vector.push_back( value );
	}
}

template< class Archive >
void SERIALIZE_NESTED_VD_VECTOR( Archive & arc, utility::vector1< utility::vector1< VD > > const & vector ) {
	arc( vector.size() );
	for ( utility::vector1< VD > item : vector ) {
		SERIALIZE_VD_VECTOR( arc, item );
	}
}

template< class Archive >
void DESERIALIZE_NESTED_VD_VECTOR( Archive & arc, utility::vector1< utility::vector1< VD > > & vector, std::map< VD, VD > const & translation ) {
	vector.clear();
	core::Size nentries; arc( nentries );
	for ( core::Size ii(1); ii <= nentries; ++ii ) {
		utility::vector1< VD > value;
		DESERIALIZE_VD_VECTOR( arc, value, translation );
		vector.push_back( value );
	}
}

template< class Archive >
void SERIALIZE_VD_VD_VECTOR_MAP( Archive & arc, std::map< VD, utility::vector1< VD > > const & map ) {
	arc( map.size() );
	for ( std::pair< VD, utility::vector1< VD > > item : map ) {
		SERIALIZE_VD( arc, item.first );
		SERIALIZE_VD_VECTOR( arc, item.second );
	}
}

template< class Archive >
void DESERIALIZE_VD_VD_VECTOR_MAP( Archive & arc, std::map< VD, utility::vector1< VD > > & map, std::map< VD, VD > const & translation ) {
	map.clear();
	core::Size nentries; arc( nentries );
	for ( core::Size ii(1); ii <= nentries; ++ii ) {
		VD key;
		utility::vector1< VD > value;
		DESERIALIZE_VD( arc, key, translation );
		DESERIALIZE_VD_VECTOR( arc, value, translation );
		map[ key ] = value;
	}
}

}
}

#endif // SERIALIZATION

#endif
