// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyz.functions.PyRosetta.hh
/// @brief  explicit instantiation of some numeric template functions for PyRosetta
/// @author Sergey Lyskov


#ifndef INCLUDED_numeric_xyz_functions_PyRosetta_HH
#define INCLUDED_numeric_xyz_functions_PyRosetta_HH

#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

namespace numeric {

#ifdef PYROSETTA
inline void _rotation_axis_angle_instantiation_for_PyRosetta()
{
	xyzMatrix_double m;
	rotation_axis_angle(m);
}
#endif

} // namespace numeric

#endif // INCLUDED_numeric_xyz_functions_PyRosetta_HH
