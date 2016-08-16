// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/PointGraphData.cc
/// @brief  Serlialization routines for the two types templated upon to define the PointGraph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/conformation/PointGraphData.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>


/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::PointGraphVertexData::save( Archive & arc ) const {
	arc( CEREAL_NVP( xyz_ ) ); // numeric::xyzVector<core::Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::PointGraphVertexData::load( Archive & arc ) {
	arc( xyz_ ); // numeric::xyzVector<core::Real>
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::PointGraphVertexData );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::PointGraphEdgeData::save( Archive & arc ) const {
	arc( CEREAL_NVP( dsq_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::PointGraphEdgeData::load( Archive & arc ) {
	arc( dsq_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::PointGraphEdgeData );
#endif // SERIALIZATION
