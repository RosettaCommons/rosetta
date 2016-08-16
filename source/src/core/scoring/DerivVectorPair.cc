// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/DerivVectorPair.hh
/// @brief  Serialization routines for the DerivVectorPair class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <core/scoring/DerivVectorPair.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::DerivVectorPair::save( Archive & arc ) const {
	arc( CEREAL_NVP( f1_ ) ); // Vector
	arc( CEREAL_NVP( f2_ ) ); // Vector
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::DerivVectorPair::load( Archive & arc ) {
	arc( f1_ ); // Vector
	arc( f2_ ); // Vector
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::DerivVectorPair );
#endif // SERIALIZATION
