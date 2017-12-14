// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/util.cc
/// @brief   Definitions for ring-related utility functions.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/chemical/rings/util.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/vector1.hh>

// Basic header
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/geometry/ring_plane.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>

// C++ header
#include <cmath>


// Construct tracer.
static basic::Tracer TR( "core.chemical.rings.util" );


namespace core {
namespace chemical {
namespace rings {

// Is the query atom axial or equatorial to the given ring or neither?
/// @details This function calculates an average plane and determines whether the coordinates of a given atom are
/// axial or equatorial to it (or neither).
/// @param   <query_atom>:      The Cartesian coordinates of the atom in question.
/// @param   <attachment_atom>: The Cartesian coordinates of the atom to which the query atom is attached, which
/// must be a member of <ring_atoms>
/// @param   <ring_atoms>:      A list of Cartesian coordinates for the points of a monocyclic ring system in
/// sequence.
/// @return  An AxEqDesignation enum type value: AXIAL, EQUATORIAL, or NEITHER
AxEqDesignation
is_atom_axial_or_equatorial_to_ring(
	Coords const & query_atom,
	Coords const & attachment_atom,
	utility::vector1< Coords > const & ring_atoms )
{
	using namespace numeric;
	using namespace numeric::constants::r;

	//JAB - Making the warnings in debug mode as we use this funciton a lot for linkage optimization.
	if ( ring_atoms.size() < 3 ) {
		TR.Debug << "A ring cannot contain fewer than 3 atoms; ";
		TR.Debug << "an axial/equatorial designation is meaningless." << std::endl;
		return NEITHER;
	}

	if ( ! ring_atoms.contains( attachment_atom ) ) {
		TR.Debug << "The attachment point for the query atom is not found in the ring; ";
		TR.Debug << "an axial/equatorial designation is meaningless." << std::endl;
		return NEITHER;
	}

	if ( ring_atoms.contains( query_atom ) ) {
		TR.Debug << "The query atom cannot be a member of the ring; ";
		TR.Debug << "an axial/equatorial designation is meaningless." << std::endl;
		return NEITHER;
	}

	// Get the vector of the bond in question.
	xyzVector< Distance > const bond( query_atom - attachment_atom );

	// Get the vector normal to the plane of best fit of the ring.
	xyzVector< Distance > const plane( geometry::vector_normal_to_ring_plane_of_best_fit( ring_atoms ) );

	// Get the angle between the two vectors.
	// Take the absolute value of the dot product, because whether the substituent is above or below the ring does
	// not matter.  This limits theta to between 0 and 90 degrees.
	Angle const theta( acos( fabs( dot( plane, bond.normalized() ) ) ) );  // in radians
	debug_assert( theta >= 0.0 && theta <= pi );

	// A "perfectly" AXIAL substituent would have a theta of 0 degrees.
	// A "perfectly" EQUATORIAL substituent would have a theta of 90 degrees, being parallel to the plane.
	// A "perfectly" NEITHER substituent would have a theta of 45 degrees.
	// An ideal chair's axial substituent would have a theta of 0 degrees,
	// and its equatorial substituent would have a theta of 180 degrees - ~110 degrees = ~70 degrees.
	// An ideal 1,4B's axial substituents on 1 and 4 would have thetas of
	// ~110 degrees - 90 degrees = ~20 degrees,
	// and its 1 and 4 equatorial substituents would have thetas of 90 degrees.
	// A planar 6-membered ring would have thetas of 90 degrees - ~110 degrees / 2 = ~35 degrees.
	// Thus, our bins need to capture the following values (in degrees):
	// AXIAL: 0, 20
	// EQUATORIAL: 70, 90
	// NEITHER: 35, 45
	// It looks like we can go with bins of:
	// AXIAL: 0 +/- 30
	// EQUATORIAL: 90 +/- 30
	// NEITHER: 45 +/- 15
	if ( theta < 30 * pi_over_180 ) {
		return AXIAL;
	} else if ( theta > 60 * pi_over_180 ) {
		// This atom could be in the middle of the ring, but that's not our problem.
		return EQUATORIAL;
	} else {
		return NEITHER;
	}
}

}  // namespace rings
}  // namespace chemical
}  // namespace core
