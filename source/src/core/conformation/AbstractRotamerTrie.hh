// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/AbstractRotamerTrie.hh
/// @brief  Base class from which all Tries will derive; this class is a handle for the tries
///         that will be defined in core/scoring, but which have to be stored by the RotamerSetBase class
///         defined in core/conformation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_conformation_AbstractRotamerTrie_hh
#define INCLUDED_core_conformation_AbstractRotamerTrie_hh

// Unit Headers
#include <core/conformation/AbstractRotamerTrie.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {

class AbstractRotamerTrie : public utility::pointer::ReferenceCount
{
public:
	~AbstractRotamerTrie() override = default;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_AbstractRotamerTrie )
#endif // SERIALIZATION


#endif
