// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/VoxelSetIterator.cxxtest.hh
/// @brief  test suite for protocols::match::VoxelSetIterator
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Utility headers
#include <utility/fixedsizearray1.hh>
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>
#include <protocols/match/VoxelSetIterator.hh>

// C++ headers
#include <string>
#include <iostream>

// AUTO-REMOVED #include <numeric/HomogeneousTransform.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map_fwd.hpp>


// --------------- Test Class --------------- //

using namespace protocols::match;

class VoxelSetIteratorTests : public CxxTest::TestSuite {

	public:

		typedef core::Vector Vector;
		typedef core::Size   Size;
		typedef core::Real   Real;
		typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
		typedef protocols::match::Real6 Real6;
	  typedef numeric::geometry::hashing::Real3 Real3;

	private:
		Vector lower;
		Vector upper;
		BoundingBox bb;

		//Real6 p;

		Real3 xyz_bin_widths;
		Real3 euler_bin_widths;
		Real3 xyz_bin_halfwidths;
		Real3 euler_bin_halfwidths;
		Size3 n_xyz_bins;
		Size3 n_euler_bins;

	public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		lower = Vector( 2, 5, 3 );
		upper = Vector( 7, 11, 6 );
		bb = BoundingBox( lower, upper );

		std::fill( xyz_bin_widths.begin(), xyz_bin_widths.end(), 1 );
		std::fill( euler_bin_widths.begin(), euler_bin_widths.end(), 10 );
		std::fill( xyz_bin_halfwidths.begin(), xyz_bin_halfwidths.end(), 0.5 );
		std::fill( euler_bin_halfwidths.begin(), euler_bin_halfwidths.end(), 5 );

		for ( Size ii = 1; ii <= 3; ++ii ) {
			n_xyz_bins[ ii ] = static_cast< Size > ((upper( ii ) - lower( ii )) / xyz_bin_widths[ ii ]);
			//if ( lower( ii ) + xyz_bin_widths[ ii ] * n_xyz_bins[ ii ] < upper( ii ) ) ++n_xyz_bins;
			n_euler_bins[ ii ] = static_cast< Size > (( ii == 3 ? 180 : 360 ) / euler_bin_widths[ ii ]);
		}

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_center_iteration() {
		/// some point in the center of all 6 dimensions; bounding box is not exceeded, and euler
		/// angles are in the center
		Real6 p;
		p[ 1 ] = 4.2; p[ 2 ] = 5.75;  p[ 3 ] = 3.8;
		p[ 4 ] = 52; p[ 5 ] = 278;  p[ 6 ] = 63;

		VoxelSetIterator voxiter( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p );

		numeric::geometry::hashing::Bin6D bin; Size pos;
		voxiter.get_bin_and_pos( bin, pos );
		//std::cout << "Vox: bin=";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << " " << bin[ ii ];
		//std::cout << " pos= " << pos << std::endl;

		TS_ASSERT( bin[ 1 ] == 2 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 5 );
		TS_ASSERT( bin[ 5 ] == 27 );
		TS_ASSERT( bin[ 6 ] == 6 );

		++voxiter;
		voxiter.get_bin_and_pos( bin, pos );

		//voxiter.get_bin_and_pos( bin, pos );
		//std::cout << "Vox: bin=";
		//for ( Size ii = 1; ii <= 6; ++ii ) std::cout << " " << bin[ ii ];
		//std::cout << " pos= " << pos << std::endl;

		TS_ASSERT( bin[ 1 ] == 2 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 5 );
		TS_ASSERT( bin[ 5 ] == 27 );
		TS_ASSERT( bin[ 6 ] == 6 );

		++voxiter;
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 2 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 5 );
		TS_ASSERT( bin[ 5 ] == 28 );
		TS_ASSERT( bin[ 6 ] == 6 );

	}

	void test_center_iterations_pos() {

		Real6 p;
		p[ 1 ] = 4.2; p[ 2 ] = 5.1;  p[ 3 ] = 3.2;
		p[ 4 ] = 52; p[ 5 ] = 271;  p[ 6 ] = 63;

		VoxelSetIterator voxiter( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p );

		for ( Size ii = 1; ii <= 64; ++ii ) {
			TS_ASSERT( ! voxiter.at_end() );

			numeric::geometry::hashing::Bin6D bin; Size pos;
			voxiter.get_bin_and_pos( bin, pos );

			TS_ASSERT( pos == ii );
			++voxiter;
		}

		TS_ASSERT( voxiter.at_end() );

	}

	void test_bounding_box_edge_z() {
		// Exceed z-range when stepping to z's upper voxel.

		Real6 p;
		p[ 1 ] = 4.2; p[ 2 ] = 5.8;  p[ 3 ] = 5.8;
		p[ 4 ] = 57; p[ 5 ] = 277;  p[ 6 ] = 67;

		VoxelSetIterator voxiter( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p );

		for ( Size ii = 1; ii <= 7; ++ii ) ++voxiter;

		numeric::geometry::hashing::Bin6D bin; Size pos;
		voxiter.get_bin_and_pos( bin, pos );


		TS_ASSERT( bin[ 1 ] == 2 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 2 );
		TS_ASSERT( bin[ 4 ] == 6 );
		TS_ASSERT( bin[ 5 ] == 28 );
		TS_ASSERT( bin[ 6 ] == 7 );

		++voxiter; //This should roll past Z's upper state, right into Y's upper state
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 2 );
		TS_ASSERT( bin[ 2 ] == 1 );
		TS_ASSERT( bin[ 3 ] == 2 );
		TS_ASSERT( bin[ 4 ] == 5 );
		TS_ASSERT( bin[ 5 ] == 27 );
		TS_ASSERT( bin[ 6 ] == 6 );


	}

	void test_bounding_box_edge_y() {
		// Exceed y-range when stepping to y's upper voxel.

		Real6 p;
		p[ 1 ] = 4.7; p[ 2 ] = 10.7;  p[ 3 ] = 3.8;
		p[ 4 ] = 57; p[ 5 ] = 277;  p[ 6 ] = 67;

		VoxelSetIterator voxiter( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p );

		for ( Size ii = 1; ii <= 15; ++ii ) ++voxiter;

		numeric::geometry::hashing::Bin6D bin; Size pos;
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 2 );
		TS_ASSERT( bin[ 2 ] == 5 );
		TS_ASSERT( bin[ 3 ] == 1 );
		TS_ASSERT( bin[ 4 ] == 6 );
		TS_ASSERT( bin[ 5 ] == 28 );
		TS_ASSERT( bin[ 6 ] == 7 );

		++voxiter; //This should roll past Y's upper state, right into X's upper state
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 3 );
		TS_ASSERT( bin[ 2 ] == 5 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 5 );
		TS_ASSERT( bin[ 5 ] == 27 );
		TS_ASSERT( bin[ 6 ] == 6 );


	}

	void test_bounding_box_edge_x() {
		// Exceed x-range when stepping to x's upper voxel.

		Real6 p;
		p[ 1 ] = 6.7; p[ 2 ] = 5.7;  p[ 3 ] = 3.8;
		p[ 4 ] = 57; p[ 5 ] = 277;  p[ 6 ] = 67;

		VoxelSetIterator voxiter( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p );

		for ( Size ii = 1; ii <= 31; ++ii ) ++voxiter;

		numeric::geometry::hashing::Bin6D bin; Size pos;
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 4 );
		TS_ASSERT( bin[ 2 ] == 1 );
		TS_ASSERT( bin[ 3 ] == 1 );
		TS_ASSERT( bin[ 4 ] == 6 );
		TS_ASSERT( bin[ 5 ] == 28 );
		TS_ASSERT( bin[ 6 ] == 7 );

		++voxiter; //This should roll past X's upper state, and end up right at the end.

		TS_ASSERT( voxiter.at_end() );
	}

	void test_wrap_phi_at_360() {
		Real6 p;
		p[ 1 ] = 3.7; p[ 2 ] = 5.7;  p[ 3 ] = 3.8;
		p[ 4 ] = 356; p[ 5 ] = 277;  p[ 6 ] = 67;

		VoxelSetIterator voxiter( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p );

		numeric::geometry::hashing::Bin6D bin; Size pos;
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 35 );
		TS_ASSERT( bin[ 5 ] == 27 );
		TS_ASSERT( bin[ 6 ] == 6 );

		for ( Size ii = 1; ii <= 3; ++ii ) ++voxiter;

		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 35 );
		TS_ASSERT( bin[ 5 ] == 28 );
		TS_ASSERT( bin[ 6 ] == 7 );

		++voxiter;
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 0 );
		TS_ASSERT( bin[ 5 ] == 27 );
		TS_ASSERT( bin[ 6 ] == 6 );

	}

	void test_wrap_psi_at_360() {
		Real6 p;
		p[ 1 ] = 3.7; p[ 2 ] = 5.7;  p[ 3 ] = 3.8;
		p[ 4 ] = 277; p[ 5 ] = 356;  p[ 6 ] = 67;

		VoxelSetIterator voxiter( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p );

		numeric::geometry::hashing::Bin6D bin; Size pos;
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 27 );
		TS_ASSERT( bin[ 5 ] == 35 );
		TS_ASSERT( bin[ 6 ] == 6 );

		++voxiter;

		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 27 );
		TS_ASSERT( bin[ 5 ] == 35 );
		TS_ASSERT( bin[ 6 ] == 7 );

		++voxiter;
		voxiter.get_bin_and_pos( bin, pos );

		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 27 );
		TS_ASSERT( bin[ 5 ] == 0 );
		TS_ASSERT( bin[ 6 ] == 6 );

	}

	void test_wrap_theta_at_0() {
		Real6 p1;
		p1[ 1 ] = 3.7; p1[ 2 ] = 5.7;  p1[ 3 ] = 3.8;
		p1[ 4 ] = 277; p1[ 5 ] =  48;  p1[ 6 ] = 3;

		Real6 p2;
		p2[ 1 ] = 3.7; p2[ 2 ] = 5.7;  p2[ 3 ] = 3.8;
		p2[ 4 ] =  97; p2[ 5 ] = 228;  p2[ 6 ] = 3;

		VoxelSetIterator voxiter1( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p1 );

		VoxelSetIterator voxiter2( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p2 );

		numeric::geometry::hashing::Bin6D bin; Size pos;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 0
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );  // already wrapped!
		TS_ASSERT( bin[ 5 ] == 22 ); // already wrapped!
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 1
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 27 );
		TS_ASSERT( bin[ 5 ] == 4 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 2
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 3
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 27 );
		TS_ASSERT( bin[ 5 ] == 5 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 4
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 5
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 28 );
		TS_ASSERT( bin[ 5 ] == 4 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 6
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 0 );

		///point 2 does not wrap!

		voxiter2.get_bin_and_pos( bin, pos );

		/// state 0
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 1
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 2
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 3
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 4
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 5
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 0 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 6
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 0 );
	}

	void test_wrap_theta_at_180() {
		Real6 p1;
		p1[ 1 ] = 3.7; p1[ 2 ] = 5.7;  p1[ 3 ] = 3.8;
		p1[ 4 ] = 277; p1[ 5 ] =  48;  p1[ 6 ] = 177;

		Real6 p2;
		p2[ 1 ] = 3.7; p2[ 2 ] = 5.7;  p2[ 3 ] = 3.8;
		p2[ 4 ] =  97; p2[ 5 ] = 228;  p2[ 6 ] = 177;

		VoxelSetIterator voxiter1( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p1 );

		VoxelSetIterator voxiter2( bb, n_xyz_bins, n_euler_bins, xyz_bin_widths,
			euler_bin_widths, xyz_bin_halfwidths, euler_bin_halfwidths, p2 );

		numeric::geometry::hashing::Bin6D bin; Size pos;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 0
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 27 );
		TS_ASSERT( bin[ 5 ] == 4 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 1
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 2
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 27 );
		TS_ASSERT( bin[ 5 ] == 5 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 3
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 4
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 28 );
		TS_ASSERT( bin[ 5 ] == 4 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 5
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter1;
		voxiter1.get_bin_and_pos( bin, pos );

		/// state 6
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 28 );
		TS_ASSERT( bin[ 5 ] == 5 );
		TS_ASSERT( bin[ 6 ] == 17 );

		///point 2 does not wrap!

		voxiter2.get_bin_and_pos( bin, pos );

		/// state 0
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 1
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 2
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 3
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 9 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 4
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 5
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 22 );
		TS_ASSERT( bin[ 6 ] == 17 );

		++voxiter2;
		voxiter2.get_bin_and_pos( bin, pos );

		/// state 6
		TS_ASSERT( bin[ 1 ] == 1 );
		TS_ASSERT( bin[ 2 ] == 0 );
		TS_ASSERT( bin[ 3 ] == 0 );
		TS_ASSERT( bin[ 4 ] == 10 );
		TS_ASSERT( bin[ 5 ] == 23 );
		TS_ASSERT( bin[ 6 ] == 17 );
	}

};
