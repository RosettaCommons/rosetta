// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <vector>

#include <numeric/coordinate_fitting/FlatLookup.hh>
#include <numeric/alignment/QCPKernel.hh>

#ifndef COORDINATEARRRAY_RMSD_FLATLOOKUP_HH_
#define COORDINATEARRRAY_RMSD_FLATLOOKUP_HH_

#include <numeric/types.hh>

namespace numeric {
namespace coordinate_fitting {

template <class Real=double>
struct CoordinateArray_RMSD_FlatLookup : public FlatLookup<Real *, numeric::Size, Real>
{
	typedef numeric::alignment::QCPKernel<Real> Kernel;

	CoordinateArray_RMSD_FlatLookup(
		Real* entry_coords,
		Real* entry_rad,
		numeric::Size n_entries,
		numeric::Size coords_per_entry
	) :
		FlatLookup<Real *, numeric::Size, Real>(),
		entry_coordinates(entry_coords),
		entry_radii(entry_rad),
		num_entries(n_entries),
		coordinates_per_entry(coords_per_entry),
		entry_size(3 * coords_per_entry),
		fragment_centers(n_entries),
		query_center()
	{
		std::vector<numeric::Size> entry_indicies(num_entries);

		for ( numeric::Size i = 0; i < num_entries; i++ ) {
			entry_indicies[i] = i;
			typename Kernel::CoordMap src_coords(&entry_coordinates[i * entry_size], 3, coordinates_per_entry);
			fragment_centers[i] = src_coords.rowwise().sum() / src_coords.cols();
		}

		this->initialize(entry_indicies.begin(), entry_indicies.end());
	}

	virtual void prepare_for_query(Real *& q)
	{
		typename Kernel::CoordMap query_coords(q, 3, coordinates_per_entry);
		query_center = query_coords.rowwise().sum() / query_coords.cols();
	}

	virtual Real entry_distance(Real *& q, numeric::Size & e)
	{
		return Kernel::calc_coordinate_rmsd(
			typename Kernel::CoordMap(q, 3, coordinates_per_entry),
			query_center,
			typename Kernel::CoordMap(&entry_coordinates[e * entry_size], 3, coordinates_per_entry),
			fragment_centers[e]
		);
	}

	virtual Real entry_radius(numeric::Size & e)
	{
		return entry_radii[e];
	}

	Real* entry_coordinates;
	Real* entry_radii;

	numeric::Size num_entries, coordinates_per_entry, entry_size;

	std::vector<typename Kernel::Point> fragment_centers;
	typename Kernel::Point query_center;
};

}
}
#endif
