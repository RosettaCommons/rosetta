
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerLinks.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <core/pack/rotamer_set/RotamerLinks.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::rotamer_set::RotamerLinks::save( Archive & arc ) const {
	arc( CEREAL_NVP( links_ ) ); // utility::vector1<utility::vector1<int> >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::rotamer_set::RotamerLinks::load( Archive & arc ) {
	arc( links_ ); // utility::vector1<utility::vector1<int> >
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::rotamer_set::RotamerLinks );
CEREAL_REGISTER_TYPE( core::pack::rotamer_set::RotamerLinks )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_rotamer_set_RotamerLinks )
#endif // SERIALIZATION
