// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#pragma once

#include <cstdlib>
#include <ndarray.h>
#include <numeric/types.hh>

namespace numeric
{

namespace alignment
{

/* Calculate coordinate array rmsd between two coordinate input arrays.
* @param  first_coordinates  Array of shape [n|1, c, 3].
* @param  second_coordinates Array of shape [n|1, c, 3].
* @param  out                Output array of shape [n].
*/
template <class Real>
void coordinate_array_rmsd(
	ndarray::Array<Real, 3, 1> const & first_coordinates,
	ndarray::Array<Real, 3, 1> const & second_coordinates,
	ndarray::Array<Real, 1> out
);

/* Align src to onto and update superimpose_coordinates with resulting superposition transform.
* @param  src_coordinates         Array of shape [n, c, 3].
* @param  onto_coordinates        Array of shape [n|1, c, 3].
* @param  superimpose_coordinates Array of shape [n, c2, 3].
* @param  out                     Output array of shape [n].
*/
template <class Real>
void coordinate_array_superimpose(
	ndarray::Array<Real, 3, 1> const & src_coordinates,
	ndarray::Array<Real, 3, 1> const & onto_coordinates,
	ndarray::Array<Real, 3, 1> const & superimpose_coordinates,
	ndarray::Array<Real, 1> out
);

/*
* Calculate broadcast coordinate rmsd between two coordinate input arrays.
* @param  first_coordinates  Array of shape [a, c, 3].
* @param  second_coordinates Array of shape [b, c, 3].
* @param  out                Output array of shape [a, b].
*/
template <class Real>
void coordinate_array_broadcast_rmsd(
	ndarray::Array<Real, 3, 1> const & first_coordinates,
	ndarray::Array<Real, 3, 1> const & second_coordinates,
	ndarray::Array<Real, 2> out
);

/*
* Calculate broadcast coordinate rmsd between two coordinate input arrays, with entries given by length and start index.
* @param   coordinates_per_entry       Number of coordinates per entry to align.
* @param  first_coordinates           Array of shape [_, 3].
* @param  first_coordinate_indicies   Array of shape [a].
* @param  second_coordinates           Array of shape [_, 3].
* @param  second_coordinate_indicies   Array of shape [b].
* @param  out                Output array of shape [a, b].
*/
template <typename Index, typename Real>
void indexed_coordinate_array_broadcast_rmsd(
	Index coordinates_per_entry,
	ndarray::Array<Real, 2, 1> const & first_coordinates,
	ndarray::Array<Index, 1> const & first_coordinate_indicies,
	ndarray::Array<Real, 2, 1> const & second_coordinates,
	ndarray::Array<Index, 1> const & second_coordinate_indicies,
	ndarray::Array<Real, 2> out
);

extern template
void coordinate_array_rmsd(
	ndarray::Array<numeric::Real, 3, 1> const & first_coordinates,
	ndarray::Array<numeric::Real, 3, 1> const & second_coordinates,
	ndarray::Array<numeric::Real, 1> out
);

extern template
void coordinate_array_superimpose(
	ndarray::Array<numeric::Real, 3, 1> const & src_coordinates,
	ndarray::Array<numeric::Real, 3, 1> const & onto_coordinates,
	ndarray::Array<numeric::Real, 3, 1> const & superimpose_coordinates,
	ndarray::Array<numeric::Real, 1> out
);

extern template
void coordinate_array_broadcast_rmsd(
	ndarray::Array<numeric::Real, 3, 1> const & first_coordinates,
	ndarray::Array<numeric::Real, 3, 1> const & second_coordinates,
	ndarray::Array<numeric::Real, 2> out
);

extern template
void indexed_coordinate_array_broadcast_rmsd(
	numeric::Size coordinates_per_entry,
	ndarray::Array<numeric::Real, 2, 1> const & first_coordinates,
	ndarray::Array<numeric::Size, 1> const & first_coordinate_indicies,
	ndarray::Array<numeric::Real, 2, 1> const & second_coordinates,
	ndarray::Array<numeric::Size, 1> const & second_coordinate_indicies,
	ndarray::Array<numeric::Real, 2> out
);

extern template
void indexed_coordinate_array_broadcast_rmsd(
	numeric::SSize coordinates_per_entry,
	ndarray::Array<numeric::Real, 2, 1> const & first_coordinates,
	ndarray::Array<numeric::SSize, 1> const & first_coordinate_indicies,
	ndarray::Array<numeric::Real, 2, 1> const & second_coordinates,
	ndarray::Array<numeric::SSize, 1> const & second_coordinate_indicies,
	ndarray::Array<numeric::Real, 2> out
);

extern template
void coordinate_array_rmsd(
	ndarray::Array<float, 3, 1> const & first_coordinates,
	ndarray::Array<float, 3, 1> const & second_coordinates,
	ndarray::Array<float, 1> out
);

extern template
void coordinate_array_superimpose(
	ndarray::Array<float, 3, 1> const & src_coordinates,
	ndarray::Array<float, 3, 1> const & onto_coordinates,
	ndarray::Array<float, 3, 1> const & superimpose_coordinates,
	ndarray::Array<float, 1> out
);

extern template
void coordinate_array_broadcast_rmsd(
	ndarray::Array<float, 3, 1> const & first_coordinates,
	ndarray::Array<float, 3, 1> const & second_coordinates,
	ndarray::Array<float, 2> out
);

extern template
void indexed_coordinate_array_broadcast_rmsd(
	numeric::Size coordinates_per_entry,
	ndarray::Array<float, 2, 1> const & first_coordinates,
	ndarray::Array<numeric::Size, 1> const & first_coordinate_indicies,
	ndarray::Array<float, 2, 1> const & second_coordinates,
	ndarray::Array<numeric::Size, 1> const & second_coordinate_indicies,
	ndarray::Array<float, 2> out
);

extern template
void indexed_coordinate_array_broadcast_rmsd(
	numeric::SSize coordinates_per_entry,
	ndarray::Array<float, 2, 1> const & first_coordinates,
	ndarray::Array<numeric::SSize, 1> const & first_coordinate_indicies,
	ndarray::Array<float, 2, 1> const & second_coordinates,
	ndarray::Array<numeric::SSize, 1> const & second_coordinate_indicies,
	ndarray::Array<float, 2> out
);

class rmsd_calc { };

}
}
