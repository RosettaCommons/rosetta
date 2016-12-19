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
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Utility headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/list.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/polymorphic.hpp>

namespace core {
namespace chemical {


/// @details First serialize a boolean representing if the restype pointer is null.
/// For non-null restypes, then serialize a boolean representing one of two paths:
/// if the restype is in a global RTS or not. If it's in a GlobalRTS, store an annotation
/// of how to pull it out of the GlobalRTS. If it's not, then just serialize the ResidueType
/// Note that Cereal will handle multiple OPs to the same object sanely.
template < class Archive >
void serialize_residue_type( Archive & arc, ResidueTypeCOP ptr )
{
	//std::cout << "Serializing RT: " << restype << std::endl;
	if ( ! ptr ) {
		bool rt_is_nonnull( false );
		arc( CEREAL_NVP( rt_is_nonnull ) );
	} else {
		bool rt_is_nonnull( true );
		TypeSetMode mode( ptr->mode() );
		ResidueTypeSetCOP global( ChemicalManager::get_instance()->residue_type_set( mode ) );
		if ( global->has( ptr ) ) {
			GlobalResidueTypeSetCOP global_rts( utility::pointer::dynamic_pointer_cast< GlobalResidueTypeSet const >( global ) );
			bool in_global_rts( true );
			arc( CEREAL_NVP( rt_is_nonnull ) );
			arc( CEREAL_NVP(  in_global_rts ) );
			arc( ::cereal::make_nvp("global_rts", global_rts->name() ) );
			arc( ::cereal::make_nvp("restype_name", ptr->name() ) );
		} else {
			bool in_global_rts( false );
			arc( CEREAL_NVP( rt_is_nonnull ) );
			arc( CEREAL_NVP( in_global_rts ) );

			// Need to do de-duplication of this particular pointer
			uint32_t id = arc.registerSharedPointer( ptr.get() );
			arc( CEREAL_NVP(id) );
			if ( id & ::cereal::detail::msb_32bit ) {
				// Hasn't been saved yet - do so.
				arc( *ptr ); // Hopefully this is the appropriate way of serializing a RT directly.
			}
		}
	}
}

INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, serialize_residue_type, ResidueTypeCOP );

/// @details See comments for serialize_residue_type
template < class Archive >
void deserialize_residue_type( Archive & arc, ResidueTypeCOP & ptr )
{
	bool rt_is_nonnull( true ); arc( rt_is_nonnull );
	if ( rt_is_nonnull ) {

		bool in_global_rts( true ); arc( in_global_rts );
		if ( in_global_rts ) {
			std::string rts_name, rt_name;
			arc( rts_name, rt_name );
			ptr = chemical::ChemicalManager::get_instance()->residue_type_set( rts_name )->name_map( rt_name ).get_self_ptr();
		} else {
			uint32_t id;
			arc( id );
			if ( id & ::cereal::detail::msb_32bit ) {
				// Hasn't been loaded yet - do so.
				ResidueTypeOP mod_restype( new ResidueType(nullptr,nullptr,nullptr,nullptr) );
				// Inform the Archive the pointer which corresponds to this id
				// (put before loading to handle circular references.)
				arc.registerSharedPointer( id, mod_restype );
				arc( *mod_restype ); // Hopefully this is the appropriate way of deserializeing a RT directly
				ptr = mod_restype;
			} else {
				ptr = utility::pointer::static_pointer_cast< ResidueType const >( arc.getSharedPointer(id) );
			}
		}
	} else {
		ptr = nullptr;
	}
}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, deserialize_residue_type, ResidueTypeCOP & );

}
}

#endif // SERIALIZATION
