// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyz.serialization.cc
/// @brief  functions for serializing xyzVector and xyzMatrix
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifdef SERIALIZATION

// Unit headers
#include <numeric/geometry/BoundingBox.srlz.hh>

// Package Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.serialization.hh>

// Utility headers
#include <utility/serialization/serialization.hh>

namespace numeric {
namespace geometry {

template < class Archive >
void
save(
	Archive & arc,
	BoundingBox< numeric::xyzVector< float > > const & xyzbb
)
{
	save_bounding_box( arc, xyzbb );
}

template < class Archive >
void
save(
	Archive & arc,
	BoundingBox< numeric::xyzVector< double > > const & xyzbb
)
{
	save_bounding_box( arc, xyzbb );
}


template < class Archive >
void
load(
	Archive & arc,
	BoundingBox< numeric::xyzVector< float > > & xyzbb
)
{
	load_bounding_box( arc, xyzbb );
}


template < class Archive >
void
load(
	Archive & arc,
	BoundingBox< numeric::xyzVector< double > > & xyzbb
)
{
	load_bounding_box( arc, xyzbb );
}


EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( BoundingBox< numeric::xyzVector< float  > > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( BoundingBox< numeric::xyzVector< double > > );


}
}

#endif
