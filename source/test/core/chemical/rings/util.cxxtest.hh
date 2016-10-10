// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/chemical/rings/util.cxxtest.hh
/// @brief   Test suite for ring-related utility functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/rings/util.hh>

// Package header
# include <core/chemical/rings/AxEqDesignation.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>

// Basic header
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/geometry/ring_plane.functions.hh>

// C++ header
#include <math.h>

static THREAD_LOCAL basic::Tracer TR( "core.chemical.rings.util.cxxtest.hh" );

class RingUtilityFunctionTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp() { core_init(); }

	// Destruction
	void tearDown()
	{}


public:  // Test //////////////////////////////////////////////////////////////
	// Confirm that axial and equatorial designations are properly assigned.
	void test_is_atom_axial_or_equatorial_to_ring()
	{
		using namespace utility;
		using namespace numeric;
		using namespace numeric::geometry;
		using namespace numeric::constants::r;
		using namespace core::chemical::rings;

		// Constants
		Real const sin60( sin( 60.0 * pi_over_180 ) );  // = square root of 3 over 2
		Real const sin45over2( sin( 45.0 * pi_over_180 ) / 2 );  // = square root of 2 over 4
		Real const sin45times2( sin( 45.0 * pi_over_180 ) * 2 );  // = square root of 2


		TR << "Testing that is_atom_axial_or_equatorial_to_ring() properly designates a ring substituent as axial, "
			"equatorial, or neither." << std::endl;

		// Definitions from IUPAC:
		// * Axial:
		//     "...[B]onds to ring atoms (and molecular entities attached to such bonds) are... axial... [if]
		//     the bonds make a relatively large... angle... with the plane containing or passing closest to
		//     a majority of the ring atoms.
		//     ...[A]xial bonds are approximately parallel to the C3 axis...."
		// * Equatorial:
		//     "...[B]onds to ring atoms (and molecular entities attached to such bonds) are... equatorial... [if]
		//     the bonds make a relatively small... angle... with the plane containing or passing closest to
		//     a majority of the ring atoms.
		//     ...[E]quatorial bonds [are] approximately parallel to two of the ring bonds."

		{
			TR << " Testing chair conformation of 6-membered ring."  << std::endl;

			// 6-membered ring in a chair conformation
			vector1< Coords > ring6( 6 );
			ring6[ 1 ] = Coords( 0.0, 1.0, sin45over2 );
			ring6[ 2 ] = Coords( sin60, 0.5, 0.0 );
			ring6[ 3 ] = Coords( sin60, -0.5, sin45over2 );
			ring6[ 4 ] = Coords( 0.0, -1.0, 0.0 );
			ring6[ 5 ] = Coords( -sin60, -0.5, sin45over2 );
			ring6[ 6 ] = Coords( -sin60, 0.5, 0.0 );

			// A "bond length" for this model system is exactly 3 over the square root of 8 or
			// 3 times the square root of 2 over 4 or 3 times sin 45 over 2.

			// The C3 axis for this model system is the z axis.

			// Perfectly axial coordinates for this chair would sit at the height of a tetrahedron resting on the xy
			// plane with an edge length of root 3, which simplifies to root 2.
			// Perfectly equatorial coordinates would sit 1 additional bond length on the y axis.
			// Realistically equatorial coordinates would sit at the the point of the tetrahedron.
			Coords const axial_atom( 0.0, 1.0, sin45times2 ),
				perfectly_equatorial_atom( 0.0, 1.0 + 3 * sin45over2, sin45over2 ),
				realistically_equatorial_atom( 0.0, 2.0, 0.0 ),
				neither_atom( 0.0, 1.5, sin45times2 / 2 ),  // just making something up for this one
				attachment_atom( 0.0, 1.0, sin45over2 );

			// The query atom is not attached to the ring.
			TS_ASSERT( ! is_atom_axial_or_equatorial_to_ring( axial_atom, ZERO_VECTOR, ring6 ) );

			// The query atom is in the ring.
			TS_ASSERT( ! is_atom_axial_or_equatorial_to_ring( attachment_atom, attachment_atom, ring6 ) );

			// The query atom is in the middle of the ring.
			//TS_ASSERT( ! is_atom_axial_or_equatorial_to_ring( ZERO_VECTOR, attachment_atom, ring6 ) );

			// The query atom is perfectly axial.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring( axial_atom, attachment_atom, ring6 ), AXIAL );

			// The query atom is perfectly equatorial by definition, though not realistic.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				perfectly_equatorial_atom, attachment_atom, ring6 ), EQUATORIAL );

			// The query atom is equatorial with realistic bond angles.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				realistically_equatorial_atom, attachment_atom, ring6 ), EQUATORIAL );

			// The query atom is neither axial nor equatorial.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				neither_atom, attachment_atom, ring6 ), NEITHER );
		}

		{
			TR << " Testing boat conformation of 6-membered ring at position 1, a clearly defined case."  << std::endl;

			// 6-membered ring in a boat conformation
			vector1< Coords > ring6( 6 );
			ring6[ 1 ] = Coords( 0.0, sin45over2, 0.5 );
			ring6[ 2 ] = Coords( sin60, 0.0, 0.0 );
			ring6[ 3 ] = Coords( sin60, -3 * sin45over2, 0 );
			ring6[ 4 ] = Coords( 0.0, -sin45times2, 0.5 );
			ring6[ 5 ] = Coords( -sin60, -3 * sin45over2, 0 );
			ring6[ 6 ] = Coords( -sin60, 0.0, 0.0 );

			// A "bond length" for this model system is exactly 3 over the square root of 8 or
			// 3 times the square root of 2 over 4 or 3 times sin 45 over 2.

			// Relative to the model system above, this boat conformation is translated along the xy plane
			// such that points 2 and 6 lie on the x axis and has been rotated along the x axis
			// such that points 2, 3, 5, and 6 all lie in the xy plane.

			// Perfectly axial coordinates would sit 1 additional bond length on the z axis.
			// Realistically axial coordinates would sit at y = 0, z = 1.5.
			// Perfectly equatorial coordinates would sit 1 additional bond length on the y axis.
			Coords const perfectly_axial_atom( 0.0, sin45over2, 0.5 + 3 * sin45over2 ),
				realistically_axial_atom( 0.0, 0.0, 1.5  ),
				equatorial_atom( 0.0, sin45times2, 0.5 ),
				attachment_atom( 0.0, sin45over2, 0.5 );

			// The query atom is perfectly axial by definition, though not realistic.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				perfectly_axial_atom, attachment_atom, ring6 ), AXIAL );

			// The query atom is axial with realistic bond angles.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				realistically_axial_atom, attachment_atom, ring6 ), AXIAL );

			// The query atom is perfectly equatorial.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				equatorial_atom, attachment_atom, ring6 ), EQUATORIAL );
		}

		{
			TR << " Testing E5 envelope conformation of 5-membered ring at position 5, a clearly defined case."  << std::endl;

			// 5-membered ring (coordinates generated by Discovery Studio)
			vector1< Coords > ring5( 5 );
			ring5[ 1 ] = Coords( -1.0227, 0.1152, -0.1809 );
			ring5[ 2 ] = Coords( -0.8198, -1.3771, 0.1630 );
			ring5[ 3 ] = Coords( 0.7211, -1.5898, 0.0512 );
			ring5[ 4 ] = Coords( 1.2678, -0.2010, -0.3471 );
			ring5[ 5 ] = Coords( 0.2733, 0.7371, 0.3289 );

			// coordinates generated by Discovery Studio
			Coords const axial_atom( 0.3466, 0.6964, 1.4156 ),
				equatorial_atom( 0.3913, 1.7679, -0.0054 ),
				attachment_atom( 0.2733, 0.7371, 0.3289 );

			// The query atom is axial.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring( axial_atom, attachment_atom, ring5 ), AXIAL );

			// The query atom is equatorial.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				equatorial_atom, attachment_atom, ring5 ), EQUATORIAL );
		}

		{
			TR << " Testing 3-membered ring, a clearly defined case."  << std::endl;

			// 3-membered ring (coordinates generated by Discovery Studio)
			vector1< Coords > ring3( 3 );
			ring3[ 1 ] = Coords( -1.114, 0.281, -0.404 );
			ring3[ 2 ] = Coords( -0.120, -0.895, -0.404 );
			ring3[ 3 ] = Coords( 0.398, 0.551, -0.296 );

			// coordinates generated by Discovery Studio
			Coords const neither_atom1( -1.76, 0.449, 0.458 ),
				neither_atom2( -1.65, 0.541, -1.316  ),
				attachment_atom( -1.114, 0.281, -0.404 );

			// The query atoms are neither axial nor equatorial.
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				neither_atom1, attachment_atom, ring3 ), NEITHER );
			TS_ASSERT_EQUALS( is_atom_axial_or_equatorial_to_ring(
				neither_atom2, attachment_atom, ring3 ), NEITHER );
		}

		{
			TR << " Testing 2-membered \"ring\"."  << std::endl;

			// 2-membered "ring"
			vector1< Coords > ring2( 2 );
			ring2[ 1 ] = Coords( 0.0 );
			ring2[ 2 ] = Coords( 1.0 );  // just a line segment

			Coords const query_atom( -0.5 ), attachment_atom( 0.0 );

			// The query atom is not attached to a ring.
			TS_ASSERT( ! is_atom_axial_or_equatorial_to_ring( query_atom, attachment_atom, ring2 ) );
		}
	}

	// Confirm that axial and equatorial designations are properly assigned.
	void test_is_atom_axial_or_equatorial()
	{
		using namespace core::chemical::rings;
		using namespace core::pose;

		// This function is not implemented at the moment?
		/*
		core_init_with_additional_options( "-include_sugars" );

		TS_TRACE( "Testing that is_atom_axial_or_equatorial() properly designates a ring substituent as axial, "
		"equatorial, or neither, including its ability to find the ring atoms and attachment point." );

		Pose glucose, serine;
		make_pose_from_saccharide_sequence( glucose, "alpha-D-Glcp", "fa_standard" );
		make_pose_from_sequence( serine, "S", "fa_standard" );

		uint const O1_index( glucose.residue( 1 ).atom_index( " O1 " ) );
		uint const O6_index( glucose.residue( 1 ).atom_index( " O6 " ) );
		uint const OG_index( serine.residue( 1 ).atom_index( " OG " ) );

		// An alpha D-sugar has an axial attachement of its O1.
		TS_ASSERT_EQUALS( is_atom_axial_or_equatorial( glucose, 1, O1_index ), AXIAL );

		// The O6 oxygen is exocyclic.
		TS_ASSERT( ! is_atom_axial_or_equatorial( glucose, 1, O6_index ) );

		// A linear residue cannot have any axial or equatorial atoms.
		TS_ASSERT( ! is_atom_axial_or_equatorial( serine, 1, OG_index ) );
		*/
		TS_ASSERT( true );
	}
};  // class RingUtilityFunctionTests
