// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/util.cc
/// @brief   Declarations for ring-related utility functions.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_rings_util_HH
#define INCLUDED_core_chemical_rings_util_HH

// Package header
#include <core/chemical/rings/AxEqDesignation.hh>

// Project header
#include <core/types.hh>

// Numeric header
#include <numeric/xyzVector.hh>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace chemical {
namespace rings {

// Type definition
typedef numeric::xyzVector< Distance > Coords;

/// @brief  Return the opposite axial/equatorial designation.
AxEqDesignation opposite_designation( AxEqDesignation designation );

/// @brief  Are the query atom coordinates axial or equatorial to the given ring or neither?
AxEqDesignation is_atom_axial_or_equatorial_to_ring(
	Coords const & query_atom,
	Coords const & attachment_atom,
	utility::vector1< Coords > const & ring_atoms );

}  // namespace rings
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_rings_util_HH
