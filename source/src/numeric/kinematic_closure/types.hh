// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/kinematic_closure/types.hh
/// @brief  Define the numeric types that are commonly used in this namespace.
/// @author Kale Kundert

#ifndef INCLUDED_numeric_kinematic_closure_types_HH
#define INCLUDED_numeric_kinematic_closure_types_HH

#include <numeric/types.hh>
#include <utility/vector1.hh>

namespace numeric {
namespace kinematic_closure {

typedef utility::vector1<numeric::Size> IndexList;
typedef utility::vector1<numeric::Real> Coordinate;
typedef utility::vector1<Coordinate> CoordinateList;

typedef utility::vector1<numeric::Real> ParameterList;
typedef utility::vector1<ParameterList> ParameterMatrix;

} // end namespace kinematic_closure
} // end namespace numeric

#endif
