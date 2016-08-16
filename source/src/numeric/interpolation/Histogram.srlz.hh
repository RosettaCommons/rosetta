// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/interpolation/Histogram.srlz.hh
/// @brief  functions for serializing numeric::interpolation::Histogram instances
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_numeric_interpolation_Histogram_SRLZ_HH
#define INCLUDED_numeric_interpolation_Histogram_SRLZ_HH

#ifdef SERIALIZATION

// Package headers
#include <numeric/interpolation/Histogram.hh>

// This file declares save/load functions for serializing and  deserializing
// Histogram objects.

namespace numeric {
namespace interpolation {

template < class Archive, class X, class Y >
void
save_histogram( Archive & arc, Histogram< X, Y > const & hist )
{
	arc( hist.densities() );
	arc( hist.first_bin() );
	arc( hist.step_size() );
	arc( hist.periodic() );
	arc( hist.bin_placement() );
	arc( hist.interpolator() );
}

template < class Archive, class X, class Y >
void
load_histogram( Archive & arc, Histogram< X, Y > & hist )
{
	arc( hist.densities() );
	arc( hist.first_bin() );
	arc( hist.step_size() );
	arc( hist.periodic() );
	arc( hist.bin_placement() );
	arc( hist.interpolator() );
	hist.set_interpolator( hist.interpolator() );
}


template < class Archive >
void
save(
	Archive & arc,
	Histogram< float, float > const & hist
);


template < class Archive >
void
save(
	Archive & arc,
	Histogram< double, double > const & hist
);

template < class Archive >
void
load(
	Archive & arc,
	Histogram< float, float > & hist
);

template < class Archive >
void
load(
	Archive & arc,
	Histogram< double, double > & hist
);

}
}

#endif // SERIALIZATION

#endif // INCLUDED_numeric_interpolation_Histogram_SRLZ_HH
