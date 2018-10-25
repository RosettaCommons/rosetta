// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#pragma once

#include "ndarray.h"
#include "ndarray/pybind11.h"

#include <numeric/alignment/rmsd_calc.hh>

namespace numeric
{
namespace alignment
{

template<typename Module> void bind_rmsd_calc(Module & m)
{
	m.def("coordinate_array_rmsd", &coordinate_array_rmsd<numeric::Real> );
	m.def("coordinate_array_superimpose", &coordinate_array_superimpose<numeric::Real> );
	m.def("coordinate_array_broadcast_rmsd", &coordinate_array_broadcast_rmsd<numeric::Real> );
	m.def("indexed_coordinate_array_broadcast_rmsd", &indexed_coordinate_array_broadcast_rmsd<numeric::Size, numeric::Real> );
	m.def("indexed_coordinate_array_broadcast_rmsd", &indexed_coordinate_array_broadcast_rmsd<numeric::SSize, numeric::Real> );
};

}
}
