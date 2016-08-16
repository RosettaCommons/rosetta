// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/ResidueTypeSet.srlz.cc
/// @brief  Serialization and deserialization routines for when working with globally-accessible ResidueTypeSets
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <core/chemical/ResidueTypeSet.srlz.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Package headers
#include <core/chemical/ChemicalManager.hh>

// Utility headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/polymorphic.hpp>

namespace core {
namespace chemical {

/// @brief "Serialize" a ResidueTypeSet by serializing its name
template < class Archive >
void serialize_residue_type_set( Archive & arc, ResidueTypeSetCOP restype_set )
{
	// We don't want to serialize the ResidueTypeSet; we just want to make it possible
	// to find the appropriate ResidueTypeSet on the node where this gets deserialized.
	std::string rts_name = restype_set->name();
	arc( rts_name );
}
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, serialize_residue_type_set, ResidueTypeSetCOP );


/// @brief "Deserialize" a ResidueTypeSet by deserializing its name and then
/// obtaining a poiner to it from the ChemicalManager
template < class Archive >
void deserialize_residue_type_set( Archive & arc, ResidueTypeSetCOP & restype_set )
{
	std::string rts_name;
	arc( rts_name );
	restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( rts_name );

}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, deserialize_residue_type_set, ResidueTypeSetCOP & );


}
}


#endif // SERIALIZATION


