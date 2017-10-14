// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/interpolation/Histogram.srlz.cc
/// @brief  functions for serializing numeric::interpolation::Histogram instances
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <numeric/interpolation/Histogram.srlz.hh>

// Utility headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// This file declares save/load functions for serializing and  deserializing
// Histogram objects.

namespace numeric {
namespace interpolation {

template < class Archive >
void
save(
	Archive & arc,
	Histogram< float, float > const & hist
)
{
	save_histogram( arc, hist );
}

template < class Archive >
void
save(
	Archive & arc,
	Histogram< double, double > const & hist
)
{
	save_histogram( arc, hist );
}

template < class Archive >
void
load(
	Archive & arc,
	Histogram< float, float > & hist
)
{
	save_histogram( arc, hist );
}

template < class Archive >
void
load(
	Archive & arc,
	Histogram< double, double > & hist
)
{
	save_histogram( arc, hist );
}

typedef Histogram< float, float > his_ff;
typedef Histogram< double, double > his_dd;

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( his_ff );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( his_dd  );


}
}

#endif // SERIALIZATION

