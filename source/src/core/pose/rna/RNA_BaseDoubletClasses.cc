// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/rna/RNA_BaseDoubletClasses.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::rna::BaseStack::save( Archive & arc ) const {
	arc( CEREAL_NVP( res1 ) ); // Size
	arc( CEREAL_NVP( res2 ) ); // Size
	arc( CEREAL_NVP( orientation ) ); // Size
	arc( CEREAL_NVP( which_side ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::rna::BaseStack::load( Archive & arc ) {
	arc( res1 ); // Size
	arc( res2 ); // Size
	arc( orientation ); // Size
	arc( which_side ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::rna::BaseStack );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::rna::BasePair::save( Archive & arc ) const {
	arc( CEREAL_NVP( res1 ) ); // Size
	arc( CEREAL_NVP( res2 ) ); // Size
	arc( CEREAL_NVP( edge1 ) ); // Size
	arc( CEREAL_NVP( edge2 ) ); // Size
	arc( CEREAL_NVP( orientation ) ); // Size
	arc( CEREAL_NVP( LW_orientation ) ); // Size
	arc( CEREAL_NVP( num_hbonds ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::rna::BasePair::load( Archive & arc ) {
	arc( res1 ); // Size
	arc( res2 ); // Size
	arc( edge1 ); // Size
	arc( edge2 ); // Size
	arc( orientation ); // Size
	arc( LW_orientation ); // Size
	arc( num_hbonds ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::rna::BasePair );

#endif // SERIALIZATION
