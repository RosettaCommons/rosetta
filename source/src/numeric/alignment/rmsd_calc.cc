// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <numeric/alignment/rmsd_calc.hh>
#include <numeric/alignment/rmsd_calc.impl.hh>
#include <numeric/types.hh>

namespace numeric
{

namespace alignment
{

template
void coordinate_array_rmsd<numeric::Real>(
	ndarray::Array<numeric::Real, 3, 1> const & first_coordinates,
	ndarray::Array<numeric::Real, 3, 1> const & second_coordinates,
	ndarray::Array<numeric::Real, 1> out
);

template
void coordinate_array_superimpose<numeric::Real>(
	ndarray::Array<numeric::Real, 3, 1> const & src_coordinates,
	ndarray::Array<numeric::Real, 3, 1> const & onto_coordinates,
	ndarray::Array<numeric::Real, 3, 1> const & superimpose_coordinates,
	ndarray::Array<numeric::Real, 1> out
);

template
void coordinate_array_broadcast_rmsd<numeric::Real>(
	ndarray::Array<numeric::Real, 3, 1> const & first_coordinates,
	ndarray::Array<numeric::Real, 3, 1> const & second_coordinates,
	ndarray::Array<numeric::Real, 2> out
);

template
void indexed_coordinate_array_broadcast_rmsd<numeric::Size, numeric::Real>(
	numeric::Size coordinates_per_entry,
	ndarray::Array<numeric::Real, 2, 1> const & first_coordinates,
	ndarray::Array<numeric::Size, 1> const & first_coordinate_indicies,
	ndarray::Array<numeric::Real, 2, 1> const & second_coordinates,
	ndarray::Array<numeric::Size, 1> const & second_coordinate_indicies,
	ndarray::Array<numeric::Real, 2> out
);

template
void indexed_coordinate_array_broadcast_rmsd<numeric::SSize, numeric::Real>(
	numeric::SSize coordinates_per_entry,
	ndarray::Array<numeric::Real, 2, 1> const & first_coordinates,
	ndarray::Array<numeric::SSize, 1> const & first_coordinate_indicies,
	ndarray::Array<numeric::Real, 2, 1> const & second_coordinates,
	ndarray::Array<numeric::SSize, 1> const & second_coordinate_indicies,
	ndarray::Array<numeric::Real, 2> out
);

template
void coordinate_array_rmsd<float>(
	ndarray::Array<float, 3, 1> const & first_coordinates,
	ndarray::Array<float, 3, 1> const & second_coordinates,
	ndarray::Array<float, 1> out
);

template
void coordinate_array_superimpose<float>(
	ndarray::Array<float, 3, 1> const & src_coordinates,
	ndarray::Array<float, 3, 1> const & onto_coordinates,
	ndarray::Array<float, 3, 1> const & superimpose_coordinates,
	ndarray::Array<float, 1> out
);

template
void coordinate_array_broadcast_rmsd<float>(
	ndarray::Array<float, 3, 1> const & first_coordinates,
	ndarray::Array<float, 3, 1> const & second_coordinates,
	ndarray::Array<float, 2> out
);

template
void indexed_coordinate_array_broadcast_rmsd<numeric::Size, float>(
	numeric::Size coordinates_per_entry,
	ndarray::Array<float, 2, 1> const & first_coordinates,
	ndarray::Array<numeric::Size, 1> const & first_coordinate_indicies,
	ndarray::Array<float, 2, 1> const & second_coordinates,
	ndarray::Array<numeric::Size, 1> const & second_coordinate_indicies,
	ndarray::Array<float, 2> out
);

template
void indexed_coordinate_array_broadcast_rmsd<numeric::SSize, float>(
	numeric::SSize coordinates_per_entry,
	ndarray::Array<float, 2, 1> const & first_coordinates,
	ndarray::Array<numeric::SSize, 1> const & first_coordinate_indicies,
	ndarray::Array<float, 2, 1> const & second_coordinates,
	ndarray::Array<numeric::SSize, 1> const & second_coordinate_indicies,
	ndarray::Array<float, 2> out
);

}
}
