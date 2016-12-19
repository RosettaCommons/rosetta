// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/ResidueTypeSet.srlz.cc
/// @brief  Serialization and deserialization routines for when working with ResidueTypeSets
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <core/chemical/ResidueTypeSet.srlz.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

// Package headers
#include <core/chemical/ChemicalManager.hh>

// Utility headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/polymorphic.hpp>

namespace core {
namespace chemical {

/// @brief "Serialize" a ResidueTypeSet, being aware of the Global versions
template < class Archive >
void serialize_residue_type_set( Archive & arc, ResidueTypeSetCOP restype_set )
{
	GlobalResidueTypeSetCOP global_rts( utility::pointer::dynamic_pointer_cast< GlobalResidueTypeSet const >( restype_set ) );
	if ( global_rts != nullptr ) {
		bool is_global( true );
		arc( CEREAL_NVP( is_global ) );
		// We don't want to serialize the ResidueTypeSet; we just want to make it possible
		// to find the appropriate ResidueTypeSet on the node where this gets deserialized.
		std::string rts_name = global_rts->name();
		arc( CEREAL_NVP( rts_name ) );
	} else {
		bool is_global( false );
		PoseResidueTypeSetCOP pose_rts( utility::pointer::dynamic_pointer_cast< PoseResidueTypeSet const >( restype_set ) );
		if ( pose_rts == nullptr ) {
			utility_exit_with_message("ERROR: Tried to serialize an unsupported ResidueTypeSet!");
		}
		arc( CEREAL_NVP( is_global ) );
		arc( CEREAL_NVP( pose_rts ) ); // Using the specific pointer subtype should get around the issue with the general base class.
	}
}
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, serialize_residue_type_set, ResidueTypeSetCOP );


/// @brief "Deserialize" a ResidueTypeSet, being aware of the Global versions
template < class Archive >
void deserialize_residue_type_set( Archive & arc, ResidueTypeSetCOP & restype_set )
{
	bool is_global( false );
	arc( is_global );
	if ( is_global ) {
		std::string rts_name;
		arc( rts_name );
		restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( rts_name );
	} else {
		PoseResidueTypeSetOP pose_rts;
		arc( pose_rts );
		restype_set = pose_rts;
	}
}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, deserialize_residue_type_set, ResidueTypeSetCOP & );


}
}


#endif // SERIALIZATION



