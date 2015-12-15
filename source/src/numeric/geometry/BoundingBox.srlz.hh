// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyz.serialization.hh
/// @brief  functions for serializing BoundingBox and xyzMatrix
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_numeric_geometry_BoundingBox_SRLZ_HH
#define INCLUDED_numeric_geometry_BoundingBox_SRLZ_HH

#ifdef SERIALIZATION

// Package headers
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/xyzVector.fwd.hh>

// This file declares save/load functions (at global scope) for serializing and
// deserializing numerix::xyzVector and numeric::xyzMatrix objects.

namespace numeric {
namespace geometry {

template < class Archive, class T >
void save_bounding_box( Archive & arc, BoundingBox< T > const & bb ) {
	arc( bb.lower(), bb.upper() );
}

template < class Archive, class T >
void load_bounding_box( Archive & arc, BoundingBox< T > & bb ) {
	T lower, upper;
	arc( lower, upper );
	bb.set_lower( lower ); bb.set_upper( upper );
}

template < class Archive, class T >
void save( Archive & arc, BoundingBox< T > const & bb ) {
	save_bounding_box( arc, bb );
}

template < class Archive, class T >
void load( Archive & arc, BoundingBox< T > & bb ) {
	load_bounding_box( arc, bb );
}

template < class Archive >
void
save(
	Archive & arch,
	BoundingBox< numeric::xyzVector< float > > const & xyzbb
);

template < class Archive >
void
save(
	Archive & arch,
	BoundingBox< numeric::xyzVector< double > > const & xyzbb
);

template < class Archive >
void
load(
	Archive & arch,
	BoundingBox< numeric::xyzVector< float > > & xyzbb
);

template < class Archive >
void
load(
	Archive & arch,
	BoundingBox< numeric::xyzVector< float > > & xyzbb
);


}
}

#endif // SERIALIZATION

#endif // INCLUDED_numeric_xyz_serialization_HH
