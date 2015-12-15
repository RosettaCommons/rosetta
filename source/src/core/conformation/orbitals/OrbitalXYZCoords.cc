// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/conformation/orbitals/OrbitalXYZCoords.cc
/// @brief Created on: Jun 30, 2011
/// @author combss

// Unit headers
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>
#endif // SERIALIZATION

namespace core {
namespace conformation {
namespace orbitals {

#ifdef     SERIALIZATION

template < class Archive >
void
OrbitalXYZCoords::save( Archive & arch ) const
{
	arch( xyz_, type_ );
}

template < class Archive >
void
OrbitalXYZCoords::load( Archive & arch )
{
	arch( xyz_, type_ );
}

SAVE_AND_LOAD_SERIALIZABLE( OrbitalXYZCoords );

#endif // SERIALIZATION

}
}
}
