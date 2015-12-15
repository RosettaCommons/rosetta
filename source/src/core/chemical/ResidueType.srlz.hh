// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/ResidueType.srlz.hh
/// @brief  Serialization and deserialization routines for when working with globally-accessible ResidueTypes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_chemical_ResidueType_SRLZ_HH
#define INCLUDED_core_chemical_ResidueType_SRLZ_HH

#ifdef SERIALIZATION

// Unit headers
#include <core/chemical/ResidueType.fwd.hh>

// C++ headers
#include <list>

namespace core {
namespace chemical {

/// @brief Depending on whether or not the input ResidueType is held in one of the global ResidueTypeSets,
/// either store the name of the ResidueTypeSet and the name of the ResidueType (for later retrieval on
/// a remote node) or serialize the ResidueType itself.
template < class Archive >
void serialize_residue_type( Archive & arc, ResidueTypeCOP restype );

/// @brief "Deserialize" a ResidueType, returning a pointer to either a globally-managed ResidueType,
/// or a genuinely deserialized ResidueType object that is not globally accessible.
template < class Archive >
void deserialize_residue_type( Archive & arc, ResidueTypeCOP & restype );

template< class Archive >
void serialize_residue_type_vector( Archive & arc, ResidueTypeCOPs const & restype );

template < class Archive >
void deserialize_residue_type_vector( Archive & arc, ResidueTypeCOPs & restype );

template< class Archive >
void serialize_residue_type_list( Archive & arc, std::list< ResidueTypeCOP > const & restype );

template < class Archive >
void deserialize_residue_type_list( Archive & arc, std::list< ResidueTypeCOP > & restype );

}
}


#endif // SERIALIZATION

#endif

