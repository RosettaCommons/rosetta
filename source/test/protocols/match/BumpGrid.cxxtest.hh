// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/BumpGrid.cxxtest.hh
/// @brief  test suite for protocols::match::SixDHasher
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/match/BumpGrid.hh>

// Utility headers
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>

// Test util headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// C++ headers
#include <string>
#include <iostream>

//Auto Headers
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map_fwd.hpp>


using namespace protocols::match;


// --------------- Test Class --------------- //

class BumpGridTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.
	core::Vector lower_corner_;
	core::Vector upper_corner_;
	numeric::geometry::BoundingBox< core::Vector > bb_;

	core::Vector lower_corner2_;
	core::Vector upper_corner2_;
	numeric::geometry::BoundingBox< core::Vector > bb2_;


	// Shared initialization goes here.
	void setUp() {
		core_init();

		lower_corner_ = core::Vector( 1.0, 1.0, 0.0 );
		upper_corner_ = core::Vector( 10.0, 11.0, 12.0 );

		bb_ = numeric::geometry::BoundingBox< core::Vector >( lower_corner_, upper_corner_ );

		lower_corner2_ = core::Vector( 4.0, -1.0, 4.0 );
		upper_corner2_ = core::Vector( 10.0, 15.0, 12.0 );

		bb2_ = numeric::geometry::BoundingBox< core::Vector >( lower_corner2_, upper_corner2_ );

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_bool3D_grid_ctor() {
		Bool3DGrid bool3d;
		bool3d.set_bounding_box( bb_ );

		TS_ASSERT( bool3d.actual_bb().lower().x() == 1.0 );
		TS_ASSERT( bool3d.actual_bb().lower().y() == 1.0 );
		TS_ASSERT( bool3d.actual_bb().lower().z() == 0.0 );

		TS_ASSERT( bool3d.actual_bb().upper().x() == 17.0 );
		TS_ASSERT( bool3d.actual_bb().upper().y() == 17.0 );
		TS_ASSERT( bool3d.actual_bb().upper().z() == 16.0 );

		//std::cout << "upperx: " <<  bool3d.actual_bb().upper().x();
		//std::cout << " uppery: " <<  bool3d.actual_bb().upper().y();
		//std::cout << " upperz: " <<  bool3d.actual_bb().upper().z() << std::endl;


	}

	// --------------- Test Cases --------------- //
	void test_bool3D_grid_ctor2() {
		Bool3DGrid bool3d;
		bool3d.set_bin_width( 0.5 );
		bool3d.set_bounding_box( bb_ );

		TS_ASSERT( bool3d.actual_bb().lower().x() == 1.0 );
		TS_ASSERT( bool3d.actual_bb().lower().y() == 1.0 );
		TS_ASSERT( bool3d.actual_bb().lower().z() == 0.0 );

		TS_ASSERT( bool3d.actual_bb().upper().x() == 13.0 );
		TS_ASSERT( bool3d.actual_bb().upper().y() == 13.0 );
		TS_ASSERT( bool3d.actual_bb().upper().z() == 12.0 );

	}


	void test_bool3D_grid_with_sphere() {
		Bool3DGrid bool3d;
		bool3d.set_bounding_box( bb_ );

		bool3d.or_by_sphere_conservative( core::Vector( 4.0, 3.0, 2.0 ), 2.0 );
		TS_ASSERT( bool3d.occupied( core::Vector( 4.25, 3.25, 2.25 ) ) );
		TS_ASSERT( bool3d.occupied( core::Vector( 3.75, 2.75, 1.75 ) ) );
		TS_ASSERT( ! bool3d.occupied( core::Vector( 5.5, 4.5, 2.5 ) ) );
		// should hash into voxel with corners 5.0 (+1), 4.0 (+1) 2.0 (+1) which would not be fully covered
		// by the sphere with radius 2
	}

	void test_bool3D_grid_with_sphere_liberal() {
		Bool3DGrid bool3d;
		bool3d.set_bounding_box( bb_ );

		core::Vector center( 4.5, 3.5, 2.0 );
		core::Real radius( 1.1 ), r2 = radius * radius;
		bool3d.or_by_sphere_liberal( center,  1.1);
		TS_ASSERT( bool3d.occupied( core::Vector( 4.25, 3.25, 2.25 ) ) );
		// The sphere touches the facet defined by the four corners:
		/// ( 4.0, 3.0, 3.0 )
		/// ( 5.0, 3.0, 3.0 )
		/// ( 4.0, 4.0, 3.0 )
		/// ( 5.0, 4.0, 3.0 )
		/// without actually touching the corners themselves.
		core::Vector corn1( 4.0, 3.0, 3.0 );
		core::Vector corn2( 5.0, 3.0, 3.0 );
		core::Vector corn3( 4.0, 4.0, 3.0 );
		core::Vector corn4( 5.0, 4.0, 3.0 );

		TS_ASSERT( center.distance_squared( corn1 ) > r2 );
		TS_ASSERT( center.distance_squared( corn2 ) > r2 );
		TS_ASSERT( center.distance_squared( corn3 ) > r2 );
		TS_ASSERT( center.distance_squared( corn4 ) > r2 );

		// In spite of the corners not being contained by the sphere, the voxel with this face
		// should still be counted as occupied.
		TS_ASSERT( bool3d.occupied( core::Vector( 4.5, 3.5, 3.5 ) ) );


		TS_ASSERT( ! bool3d.occupied( core::Vector( 4.5, 3.5, 4.5 ) ) );
		// should hash into voxel with corners 5.0 (+1), 4.0 (+1) 2.0 (+1) which would not be fully covered
		// by the sphere with radius 2
	}

	void test_bool3D_grid_with_sphere2() {
		Bool3DGrid bool3d1;
		bool3d1.set_bounding_box( bb_ );

		bool3d1.or_by_sphere_conservative( core::Vector( 4.0, 2.5, 3.5 ), 0.75 );
		bool3d1.or_by_sphere_conservative( core::Vector( 5.0, 2.5, 3.5 ), 0.75 );

		TS_ASSERT( ! bool3d1.occupied( core::Vector( 4.5, 2.5, 3.5 ) ));

		Bool3DGrid bool3d2;
		bool3d2.set_bounding_box( bb_ );

		utility::vector1< std::pair< core::Vector, core::Real > > spheres; spheres.reserve( 2 );
		spheres.push_back( std::make_pair( core::Vector( 4.0, 2.5, 3.5 ), 0.75*0.75 ));
		spheres.push_back( std::make_pair( core::Vector( 5.0, 2.5, 3.5 ), 0.75*0.75 ));

		bool3d2.or_by_spheres_conservative( spheres );

		TS_ASSERT( bool3d2.occupied( core::Vector( 4.5, 2.5, 3.5 ) ));
	}

	void test_bool3D_create_shifted_grid() {
		Bool3DGrid bool3d1;
		bool3d1.set_bounding_box( bb_ );

		core::Vector new_low_corner( 4.5, 2.5, 3.5 );
		core::Vector new_upper_corner( 8.5, 4.5, 5.5 );

		Bool3DGrid bool3d2 = bool3d1.create_grid_for_bb(
			numeric::geometry::BoundingBox< core::Vector >( new_low_corner, new_upper_corner ));
		TS_ASSERT( bool3d2.actual_bb().lower().x() == 4.0 );
		TS_ASSERT( bool3d2.actual_bb().lower().y() == 2.0 );
		TS_ASSERT( bool3d2.actual_bb().lower().z() == 3.0 );

		TS_ASSERT( bool3d2.actual_bb().upper().x() == 12.0 );
		TS_ASSERT( bool3d2.actual_bb().upper().y() == 10.0 );
		TS_ASSERT( bool3d2.actual_bb().upper().z() == 11.0 );
	}

	void test_bool3D_corners() {
		Bool3DGrid bool3d;
		bool3d.set_bounding_box( bb_ );

		numeric::geometry::hashing::Bin3D bin;
		bin[ 1 ] = 2;
		bin[ 2 ] = 2;
		bin[ 3 ] = 3;

		Bool3DGrid::CornerPoints corners1 = bool3d.corners( bin );
		TS_ASSERT( corners1[ 1 ].distance_squared( core::Vector( 3, 3, 3 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 2 ].distance_squared( core::Vector( 3, 3, 4 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 3 ].distance_squared( core::Vector( 3, 4, 3 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 4 ].distance_squared( core::Vector( 3, 4, 4 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 5 ].distance_squared( core::Vector( 4, 3, 3 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 6 ].distance_squared( core::Vector( 4, 3, 4 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 7 ].distance_squared( core::Vector( 4, 4, 3 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 8 ].distance_squared( core::Vector( 4, 4, 4 ) ) < 1e-6 );
	}

	void test_bool3D_corners2() {
		Bool3DGrid bool3d;
		bool3d.set_bin_width( 0.5 );
		bool3d.set_bounding_box( bb_ );

		numeric::geometry::hashing::Bin3D bin;
		bin[ 1 ] = 4;
		bin[ 2 ] = 4;
		bin[ 3 ] = 6;

		Bool3DGrid::CornerPoints corners1 = bool3d.corners( bin );
		TS_ASSERT( corners1[ 1 ].distance_squared( core::Vector( 3,   3,   3   ) ) < 1e-6 );
		TS_ASSERT( corners1[ 2 ].distance_squared( core::Vector( 3,   3,   3.5 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 3 ].distance_squared( core::Vector( 3,   3.5, 3   ) ) < 1e-6 );
		TS_ASSERT( corners1[ 4 ].distance_squared( core::Vector( 3,   3.5, 3.5 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 5 ].distance_squared( core::Vector( 3.5, 3,   3   ) ) < 1e-6 );
		TS_ASSERT( corners1[ 6 ].distance_squared( core::Vector( 3.5, 3,   3.5 ) ) < 1e-6 );
		TS_ASSERT( corners1[ 7 ].distance_squared( core::Vector( 3.5, 3.5, 3   ) ) < 1e-6 );
		TS_ASSERT( corners1[ 8 ].distance_squared( core::Vector( 3.5, 3.5, 3.5 ) ) < 1e-6 );
	}

	void test_bool3D_or_overlap() {
		Bool3DGrid bool3d1, bool3d2;
		bool3d1.set_bounding_box( bb_ );
		bool3d2.set_bounding_box( bb2_ );

		TS_ASSERT( bool3d2.actual_bb().lower().x() ==  4.0 );
		TS_ASSERT( bool3d2.actual_bb().lower().y() == -1.0 );
		TS_ASSERT( bool3d2.actual_bb().lower().z() ==  4.0 );

		TS_ASSERT( bool3d2.actual_bb().upper().x() == 12.0 );
		TS_ASSERT( bool3d2.actual_bb().upper().y() == 15.0 );
		TS_ASSERT( bool3d2.actual_bb().upper().z() == 12.0 );


		core::Vector center( 7.5, 5.5, 8.5 );
		core::Vector contained( 8.5, 5.5, 8.5 );
		core::Vector edge_case( 9.5, 5.5, 8.5 );

		bool3d2.or_by_sphere_conservative( center, 2 );

		TS_ASSERT( bool3d2.occupied( center ) );
		TS_ASSERT( bool3d2.occupied( contained ) );
		TS_ASSERT( !bool3d2.occupied( edge_case ) );
		TS_ASSERT( !bool3d1.occupied( center ) );
		TS_ASSERT( !bool3d1.occupied( contained ) );
		TS_ASSERT( !bool3d1.occupied( edge_case ) );

		bool3d1.or_with( bool3d2 );
		TS_ASSERT( bool3d1.occupied( center ) );
		TS_ASSERT( bool3d1.occupied( contained ) );
		TS_ASSERT( !bool3d1.occupied( edge_case ) );

	}

	void test_bool3D_sphere_coverage() {
		typedef core::Size Size;

		Bool3DGrid bool3d1;
		bool3d1.set_bin_width( 0.25 );
		bool3d1.set_bounding_box( bb_ );

		numeric::geometry::hashing::Bin3D ndims = bool3d1.dimsizes();
		numeric::geometry::hashing::Bin3D bin;

		/// Grid should start out empty
		for ( Size ii = 0; ii < ndims[ 1 ]; ++ii ) {
			bin[ 1 ] = ii;
			for ( Size jj = 0; jj < ndims[ 2 ]; ++jj ) {
				bin[ 2 ] = jj;
				for ( Size kk = 0; kk < ndims[ 3 ]; ++kk ) {
					bin[ 3 ] = kk;
					Bool3DGrid::CornerPoints c = bool3d1.corners( bin );
					TS_ASSERT( !bool3d1.occupied( bin ) );
				}
			}
		}


		core::Vector const center( 4.0, 2.5, 3.5 );
		core::Real const radius( 1.75 );
		core::Real const rad2( radius * radius );

		bool3d1.or_by_sphere_conservative( center, radius );

		//Bool3DGridKinemageWriter writer;
		//writer.set_write_facets( true );
		//writer.write_grid_to_kinemage( std::cout, "test1", bool3d1 );

		for ( Size ii = 0; ii < ndims[ 1 ]; ++ii ) {
			bin[ 1 ] = ii;
			for ( Size jj = 0; jj < ndims[ 2 ]; ++jj ) {
				bin[ 2 ] = jj;
				for ( Size kk = 0; kk < ndims[ 3 ]; ++kk ) {
					bin[ 3 ] = kk;
					Bool3DGrid::CornerPoints c = bool3d1.corners( bin );
					if ( bool3d1.occupied( bin ) ) {
						bool all_in_sphere( true );
						for ( Size ll = 1; ll <= 8; ++ll ) {
							if ( c[ ll ].distance_squared( center ) > rad2 ) {
								std::cerr << "Corner not contained within sphere! " << c[ ll ].x() << " " << c[ ll ].y() << " " << c[ ll ].z() << " bin: " << bin[ 1 ] << " " << bin[ 2 ] << " " << bin[ 3 ] << " " << c[ ll ].distance( center ) << std::endl;
								core::Vector voxel_center = bool3d1.bin_center( bin );
								std::cerr << "Center at " << voxel_center.x() << " " << voxel_center.y() << " " << voxel_center.z() << " dist: " << voxel_center.distance( center ) << " ?< " << 1.75 - 0.125 * std::sqrt( 3 ) << std::endl;
								all_in_sphere = false;
							}
						}
						TS_ASSERT( all_in_sphere );
					} else {
						bool one_outside_sphere( false );
						for ( Size ll = 1; ll <= 8; ++ll ) {
							if ( c[ ll ].distance_squared( center ) > rad2 ) {
								one_outside_sphere = true;
							}
						}
						if ( ! one_outside_sphere ) {
							std::cerr << "All corners contained within sphere?! " << c[ 1 ].x() << " " << c[ 1 ].y() << " " << c[ 1 ].z() << " bin: " << bin[ 1 ] << " " << bin[ 2 ] << " " << bin[ 3 ] << std::endl;
						}
						TS_ASSERT( one_outside_sphere );
					}

				}
			}
		}

	}

	void test_bool3D_write_kinemage() {
		typedef core::Size Size;

		Bool3DGrid bool3d1;
		bool3d1.set_bin_width( 1.0 );
		bool3d1.set_bounding_box( bb_ );

		numeric::geometry::hashing::Bin3D ndims = bool3d1.dimsizes();
		numeric::geometry::hashing::Bin3D bin;

		/// Grid should start out empty
		for ( Size ii = 0; ii < ndims[ 1 ]; ++ii ) {
			bin[ 1 ] = ii;
			for ( Size jj = 0; jj < ndims[ 2 ]; ++jj ) {
				bin[ 2 ] = jj;
				for ( Size kk = 0; kk < ndims[ 3 ]; ++kk ) {
					bin[ 3 ] = kk;
					Bool3DGrid::CornerPoints c = bool3d1.corners( bin );
					TS_ASSERT( !bool3d1.occupied( bin ) );
				}
			}
		}


		core::Vector const center( 4.0, 2.5, 3.5 );
		core::Real const radius( 1.75 );
		//core::Real const rad2( radius * radius );

		bool3d1.or_by_sphere_conservative( center, radius );

		Bool3DGridKinemageWriter writer;
		writer.set_write_facets( true );
		//writer.write_grid_to_kinemage( std::cout, "test1", bool3d1 );

		std::ostringstream sstream;
		writer.write_grid_to_kinemage( sstream, "test1", bool3d1 );
		//std::cout << sstream.str() << std::endl;

		std::string correct_kinemage =
			"@group {test1} dominant\n"
			"@vectorlist {} color= red width= 1\n"
			"{\"}P U 3.75 2.75 3.75\n"
			"{\"} U 3.25 2.75 3.75\n"
			"{\"} U 3.25 2.25 3.75\n"
			"{\"} U 3.75 2.25 3.75\n"
			"{\"} U 3.75 2.75 3.75\n"
			"{\"} U 3.75 2.75 3.25\n"
			"{\"} U 3.25 2.75 3.25\n"
			"{\"} U 3.25 2.25 3.25\n"
			"{\"} U 3.75 2.25 3.25\n"
			"{\"} U 3.75 2.75 3.25\n"
			"{\"}P U 3.25 2.75 3.75\n"
			"{\"} U 3.25 2.75 3.25\n"
			"{\"}P U 3.25 2.25 3.75\n"
			"{\"} U 3.25 2.25 3.25\n"
			"{\"}P U 3.75 2.25 3.75\n"
			"{\"} U 3.75 2.25 3.25\n"
			"{\"}P U 4.75 2.75 3.75\n"
			"{\"} U 4.25 2.75 3.75\n"
			"{\"} U 4.25 2.25 3.75\n"
			"{\"} U 4.75 2.25 3.75\n"
			"{\"} U 4.75 2.75 3.75\n"
			"{\"} U 4.75 2.75 3.25\n"
			"{\"} U 4.25 2.75 3.25\n"
			"{\"} U 4.25 2.25 3.25\n"
			"{\"} U 4.75 2.25 3.25\n"
			"{\"} U 4.75 2.75 3.25\n"
			"{\"}P U 4.25 2.75 3.75\n"
			"{\"} U 4.25 2.75 3.25\n"
			"{\"}P U 4.25 2.25 3.75\n"
			"{\"} U 4.25 2.25 3.25\n"
			"{\"}P U 4.75 2.25 3.75\n"
			"{\"} U 4.75 2.25 3.25\n"
			"@ribbonlist {} color= pinktint master= {grid facets} \n"
			"{\"} U 3.75 2.25 3.25\n"
			"{\"} U 3.75 2.25 3.75\n"
			"{\"} U 3.75 2.75 3.25\n"
			"{\"} U 3.75 2.75 3.75\n"
			"{\"} U 3.25 2.75 3.25\n"
			"{\"} U 3.25 2.75 3.75\n"
			"{\"} U 3.25 2.25 3.25\n"
			"{\"} U 3.25 2.25 3.75\n"
			"{\"} U 3.75 2.25 3.25\n"
			"{\"} U 3.75 2.25 3.75\n"
			"@ribonlist {} color= pinktint master= {grid facets} \n"
			"{\"} U 3.75 2.75 3.25\n"
			"{\"}  U 3.75 2.25 3.25\n"
			"{\"} U 3.25 2.75 3.25\n"
			"{\"} U 3.25 2.25 3.25\n"
			"@ribonlist {} color= pinktint master= {grid facets} \n"
			"{\"} U 3.75 2.75 3.75\n"
			"{\"} U 3.25 2.75 3.75\n"
			"{\"} U 3.75 2.25 3.75\n"
			"{\"} U 3.25 2.25 3.75\n"
			"@ribbonlist {} color= pinktint master= {grid facets} \n"
			"{\"} U 4.75 2.25 3.25\n"
			"{\"} U 4.75 2.25 3.75\n"
			"{\"} U 4.75 2.75 3.25\n"
			"{\"} U 4.75 2.75 3.75\n"
			"{\"} U 4.25 2.75 3.25\n"
			"{\"} U 4.25 2.75 3.75\n"
			"{\"} U 4.25 2.25 3.25\n"
			"{\"} U 4.25 2.25 3.75\n"
			"{\"} U 4.75 2.25 3.25\n"
			"{\"} U 4.75 2.25 3.75\n"
			"@ribonlist {} color= pinktint master= {grid facets} \n"
			"{\"} U 4.75 2.75 3.25\n"
			"{\"}  U 4.75 2.25 3.25\n"
			"{\"} U 4.25 2.75 3.25\n"
			"{\"} U 4.25 2.25 3.25\n"
			"@ribonlist {} color= pinktint master= {grid facets} \n"
			"{\"} U 4.75 2.75 3.75\n"
			"{\"} U 4.25 2.75 3.75\n"
			"{\"} U 4.75 2.25 3.75\n"
			"{\"} U 4.25 2.25 3.75\n";
		TS_ASSERT( sstream.str() == correct_kinemage );
	}

	void test_bumpgrid_for_pose() {
		using namespace core;

		core::pose::Pose trpcage = create_trpcage_ideal_pose();

		//clock_t starttime = clock();

		BumpGridOP bgop = bump_grid_to_enclose_pose( trpcage );

		Bool3DGridKinemageWriter writer;
		writer.set_skip_completely_buried_positions( true );
		//writer.write_grid_to_kinemage( std::cout, "test1", bool3d1 );

		//for ( Size loop = 1; loop <= 100; ++loop ) {

		for ( Size ii = 1; ii <= trpcage.size(); ++ii ) {
			BumpGridOP resbgop = bump_grid_to_enclose_residue_backbone( trpcage.residue( ii ), *bgop );
			fill_grid_with_backbone_heavyatom_spheres( trpcage.residue( ii ), *resbgop );
			//std::string fname = "trpcage_bumpgrid_res_" + utility::to_string( ii ) + ".kin";
			//std::ofstream outfile;
			//outfile.open( fname.c_str() );
			//outfile << "kinemage 1\n";
			//writer.write_grid_to_kinemage( outfile, "C_ALA", resbgop->grid( C_ALA ) );
			//writer.write_grid_to_kinemage( outfile, "OXY", resbgop->grid( OXY ) );
			bgop->or_with( *resbgop );
		}

		//}

		//clock_t stoptime = clock();
		//std::cout << "Constructed bbgrid on average in " << ((double) stoptime - starttime ) / ( 100 * CLOCKS_PER_SEC ) << " second" << std::endl;


		//std::string fname = "trpcage_bumpgrid.kin";
		//std::ofstream outfile;
		//outfile.open( fname.c_str() );
		//outfile << "kinemage 1\n";
		//writer.write_grid_to_kinemage( outfile, "C_ALA", bgop->grid( C_ALA ) );
		//writer.write_grid_to_kinemage( outfile, "OXY", bgop->grid( OXY ) );

	}

};
