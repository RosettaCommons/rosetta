// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/ResidueType.srlz.cc
/// @brief  Serialization and deserialization routines for when working with globally-accessible ResidueTypes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <core/chemical/ResidueType.srlz.hh>

// Package headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Utility headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/list.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/polymorphic.hpp>

namespace core {
namespace chemical {


/// @details First serialize a boolean representing which of the two paths
/// this code will go down, then serialize the data for those two paths.
/// @note The current signal for whether a ResidueType is a member of a
/// globally-held ResidueTypeSet is whether or not the ResidueType is held
/// in a ResidueTypeSet at all.  This will fail if a non-global ResidueTypeSet
/// is used to manage a ResidueType.
template < class Archive >
void serialize_residue_type( Archive & arc, ResidueTypeCOP restype )
{
	//std::cout << "Serializing RT: " << restype << std::endl;
	if ( ! restype ) {
		bool rt_is_nonnull( false );
		arc( rt_is_nonnull );
	} else if ( restype->in_residue_type_set() ) {
		bool rt_is_nonnull( true );
		bool yes_restype_is_in_rts( true );
		arc( rt_is_nonnull );
		arc( yes_restype_is_in_rts );
		arc( restype->residue_type_set()->name() );
		arc( restype->name() );
	} else {
		// TEMPORARY! arc( restype );
		std::cout << "UNIMPLEMENTED CASE IN ResidueType.srlz.cc::serialize_residue_type!" << std::endl;
	}
}
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, serialize_residue_type, ResidueTypeCOP );

/// @details See comments for serialize_residue_type
template < class Archive >
void deserialize_residue_type( Archive & arc, ResidueTypeCOP & restype )
{
	bool rt_is_nonnull( true ); arc( rt_is_nonnull );
	if ( rt_is_nonnull ) {

		bool in_global_rts( true ); arc( in_global_rts );
		if ( in_global_rts ) {
			std::string rts_name, rt_name;
			arc( rts_name, rt_name );
			restype = chemical::ChemicalManager::get_instance()->residue_type_set( rts_name )->name_map( rt_name ).get_self_ptr();
		} else {
			// TEMP!
			std::cout << "UNIMPLEMENTED CASE IN ResidueType.srlz.cc::deserialize_residue_type!" << std::endl;
		}
	} else {
		restype = 0;
	}
}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, deserialize_residue_type, ResidueTypeCOP & );

template< class Archive >
void serialize_residue_type_vector( Archive & arc, ResidueTypeCOPs const & restypes )
{
	arc( restypes.size() );
	for ( Size ii = 1; ii <= restypes.size(); ++ii ) {
		serialize_residue_type( arc, restypes[ ii ] );
	}
}
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, serialize_residue_type_vector, ResidueTypeCOPs const & );


template < class Archive >
void deserialize_residue_type_vector( Archive & arc, ResidueTypeCOPs & restypes )
{
	Size size_of_array( 0 ); arc( size_of_array );
	restypes.resize( size_of_array );
	for ( Size ii = 1; ii <= size_of_array; ++ii ) {
		deserialize_residue_type( arc, restypes[ ii ] );
	}
}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, deserialize_residue_type_vector, ResidueTypeCOPs & );

template< class Archive >
void serialize_residue_type_list( Archive & arc, std::list< ResidueTypeCOP > const & restypes )
{
	arc( restypes.size() );
	for ( auto iter = restypes.begin(), iter_end = restypes.end(); iter != iter_end; ++iter ) {
		serialize_residue_type( arc, *iter );
	}
}
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, serialize_residue_type_list, std::list< ResidueTypeCOP > const &  );


template < class Archive >
void deserialize_residue_type_list( Archive & arc, std::list< ResidueTypeCOP > & restypes )
{
	Size numrestypes; arc( numrestypes );
	restypes.clear();
	for ( Size ii = 1; ii <= numrestypes; ++ii ) {
		ResidueTypeCOP iirt;
		deserialize_residue_type( arc, iirt );
		restypes.push_back( iirt );
	}

}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, deserialize_residue_type_list, std::list< ResidueTypeCOP > & );


}
}

#endif // SERIALIZATION
