// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Header file for shared types.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_kinematic_closure_types_HH
#define INCLUDED_protocols_kinematic_closure_types_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <numeric/conversions.hh>
#include <numeric/kinematic_closure/types.hh>
#include <numeric/kinematic_closure/vector.hh>

namespace protocols {
namespace kinematic_closure {

using core::Size;
using core::Real;
using core::pose::Pose;

using protocols::loops::Loop;

using numeric::conversions::AngleUnit;
using numeric::kinematic_closure::IndexList;
using numeric::kinematic_closure::Coordinate;
using numeric::kinematic_closure::CoordinateList;
using numeric::kinematic_closure::ParameterList;
using numeric::kinematic_closure::ParameterMatrix;
using numeric::kinematic_closure::operator <<;

} // end namespace kinematic_closure
} // end namespace protocols

#endif
