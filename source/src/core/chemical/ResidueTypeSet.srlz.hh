// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/ResidueTypeSet.srlz.hh
/// @brief  Serialization and deserialization routines for when working with globally-accessible ResidueTypeSets
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_chemical_ResidueTypeSet_SRLZ_HH
#define INCLUDED_core_chemical_ResidueTypeSet_SRLZ_HH

#ifdef SERIALIZATION

// Unit headers
#include <core/chemical/ResidueTypeSet.fwd.hh>

namespace core {
namespace chemical {

/// @brief "Serialize" a ResidueTypeSet by serializing its name
template < class Archive >
void serialize_residue_type_set( Archive & arc, ResidueTypeSetCOP restype );

/// @brief "Deserialize" a ResidueTypeSet by deserializing its name and then
/// obtaining a poiner to it from the ChemicalManager
template < class Archive >
void deserialize_residue_type_set( Archive & arc, ResidueTypeSetCOP & restype );


}
}


#endif // SERIALIZATION

#endif

