// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/sasa.cxxtest.hh
/// @brief  test suite for core::scoring::sasa
/// @author Ron Jacak
/// @author Jared Adolf-Bryfogle

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/sasa.hh>
#include <core/scoring/sasa/util.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
 // for core::pose::make_pose_from_sequence

// AUTO-REMOVED #include <core/scoring/types.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <utility/vector1.hh>



// --------------- Test Class --------------- //

// using declarations
using namespace core;

class SasaTests : public CxxTest::TestSuite {

	public:

	Real TOLERATED_ERROR;
	pose::Pose pose;
	scoring::ScoreFunctionOP sf;

	Real total_sasa;
	id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;
	Real probe_radius;

	id::AtomID_Map< bool > atom_subset;

	int phi_index, theta_index;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();
		TOLERATED_ERROR = 0.1;

		// read sasa datafiles
		core::scoring::input_sasa_dats();

		total_sasa = 0.0;
		atom_sasa.clear();
		rsd_sasa.clear();
		probe_radius = 1.4;
		atom_subset.clear();

	}


	// Shared finalization goes here.
	void tearDown() {}


	// --------------- Test Cases --------------- //
	void test_get_overlap() {

		int degree_of_overlap;

		// void get_overlap( Real const radius_a, Real const radius_b, Real const distance_ijxyz, int & degree_of_overlap );
		// tests the get_overlap() method of sasa.cc.  The specific cases we want to test include the following:
		// 1) distance between 'a' and 'b' less than an internal function variable epsilon
		//    a) radius a < radius b, will return 100
		//    b) radius b < radius a, 1
		core::scoring::sasa::get_legrand_atomic_overlap( 1.40, 2.20, 0.005, degree_of_overlap );
		TS_ASSERT_DELTA( degree_of_overlap, 100, TOLERATED_ERROR );

		core::scoring::sasa::get_legrand_atomic_overlap( 2.20, 1.40, 0.005, degree_of_overlap );
		TS_ASSERT_DELTA( degree_of_overlap, 1, TOLERATED_ERROR );


		// 2) radius_b + distance_ijxyz <= radius_a, means b completely engulfed by a so return overlap on a of 1
		core::scoring::sasa::get_legrand_atomic_overlap( 3.20, 2.20, 1.00, degree_of_overlap );
		TS_ASSERT_DELTA( degree_of_overlap, 1, TOLERATED_ERROR );


		// 3) radius_a + distance_ijxyz <= radius_b, means a completely engulfed by b so return overlap on a of 100
		core::scoring::sasa::get_legrand_atomic_overlap( 2.20, 3.20, 1.00, degree_of_overlap );
		TS_ASSERT_DELTA( degree_of_overlap, 100, TOLERATED_ERROR );


		// 4) some cases not in the three above
		//    a) make sure to have one where theta is very small < 1degree
		//    b) theta very big (close to 180deg)
		core::scoring::sasa::get_legrand_atomic_overlap( 2.20, 2.20, 0.15, degree_of_overlap );
		TS_ASSERT_DELTA( degree_of_overlap, 49, TOLERATED_ERROR );

		core::scoring::sasa::get_legrand_atomic_overlap( 1.40, 2.20, 1.0, degree_of_overlap ); // cos theta = -0.671, 1.671*50=84
		TS_ASSERT_DELTA( degree_of_overlap, 84, TOLERATED_ERROR );

		core::scoring::sasa::get_legrand_atomic_overlap( 1.40, 2.20, 2.0, degree_of_overlap ); // cos theta = 0.2, 41
		TS_ASSERT_DELTA( degree_of_overlap, 41, TOLERATED_ERROR );

		core::scoring::sasa::get_legrand_atomic_overlap( 1.40, 2.20, 3.0, degree_of_overlap ); // my calcs give cos theta = 0.728, 14
		TS_ASSERT_DELTA( degree_of_overlap, 14, TOLERATED_ERROR );

	}


	// ronj leaving this test commented out because I could never figure out what the true values were supposed to be
	// ronj based on the output of the function and the data in the database files
	void x_test_get_orientation() {

		// lines 1-5 of the database file sampling/SASA-angles.dat are the following:
		// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
		// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
		// 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 17 17 17 17 17 17 17 17 17 17 17 17 17 27 27 27 27 27 27 27 27 27 27 27 27 37 37 37 37 37 37 37 37 37 37 37 37 37 6 6 6 6 6 6
		// 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 17 17 17 17 17 17 17 17 17 17 17 17 28 27 27 27 27 27 27 27 27 27 27 27 27 37 37 37 37 37 37 37 37 37 37 37 37 37 6 6 6 6 6 6
		// 6 6 6 7 7 7 7 7 7 7 7 5 5 5 5 5 18 18 18 18 18 18 18 17 17 17 17 17 17 28 28 28 28 28 28 28 27 27 27 27 27 27 38 38 38 38 38 38 38 37 37 37 37 37 45 45 45 45 45 45 45 45 6 6

		// so line 5, or phi bin 5 (22.5-28.125 degrees, 0.392-0.4908 rad) and theta bin 2 (5.625-11.25 degrees, 0.098-0.196 radians) give a dot location of 6
		// use calc3d to get cartesian for ( 1.0, theta: 0.14, phi: 0.42 )
		Real dist = 1.0;
		Vector a( 0.1274, 0.0569, 0.9902 );
		Vector b( 0, 0, 0 );
		core::scoring::sasa::get_legrand( a, b, phi_index, theta_index, dist );
		TS_ASSERT_DELTA( phi_index, 5, TOLERATED_ERROR );
		TS_ASSERT_DELTA( theta_index, 2, TOLERATED_ERROR );

		// according to sphere.txt, dot 6 is located at either 1) theta 14.546 deg (0.254 rad) and phi 72 deg (1.257 rad)
		// or 2) theta 165.6 deg (2.89 rad) and phi 108.3 deg (1.89 rad). neither of these are dot 6 if converted to bins and looked up in table.
		// dot_coords[6] = -1* Vector( 0.0776089, 0.238856, 0.967949 );


		// so line 5, or theta bin 5 (11.25-14.0625 degrees) and phi bin 2 (5.625-11.25 degrees, 0.098-0.196 radians) gives a dot location of 6
		// use calc3d to get cartesian for ( 1.0, theta: 0.14, phi: 0.42 )


		Vector aa( 0.4038, 0.05690, 0.9131 ); // reverse theta and phi (1.0, theta: 0.42, phi: 0.14)
		core::scoring::sasa::get_legrand_orientation( aa, b, phi_index, theta_index, dist );
		TS_ASSERT_DELTA( phi_index, 5, TOLERATED_ERROR );
		TS_ASSERT_DELTA( theta_index, 2, TOLERATED_ERROR );

		//std::cout << "\n" << std::endl;


		// lines 1-5 of the database file sampling/SASA-angles.dat are the following:
		// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
		// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
		// 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 17 17 17 17 17 17 17 17 17 17 17 17 17 27 27 27 27 27 27 27 27 27 27 27 27 37 37 37 37 37 37 37 37 37 37 37 37 37 6 6 6 6 6 6
		// 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 17 17 17 17 17 17 17 17 17 17 17 17 28 27 27 27 27 27 27 27 27 27 27 27 27 37 37 37 37 37 37 37 37 37 37 37 37 37 6 6 6 6 6 6
		// 6 6 6 7 7 7 7 7 7 7 7 5 5 5 5 5 18 18 18 18 18 18 18 17 17 17 17 17 17 28 28 28 28 28 28 28 27 27 27 27 27 27 38 38 38 38 38 38 38 37 37 37 37 37 45 45 45 45 45 45 45 45 6 6
		//
		// phi bin 4 (16.8-22.5 degrees, 0.2945-0.392rad) and theta bin 33 (180-185.6 degrees, 3.14-3.239rad) give a dot location of 28 (27 all around it).
		// let's use calc3d to give a cartesian for ( 1.0, phi: 0.350, theta: 3.19):

		dist = 1.0;
		Vector c( -0.04545, -0.01659, -0.99883 );
		Vector d( 0, 0, 0 );
		core::scoring::sasa::get_legrand_orientation( c, d, phi_index, theta_index, dist );
		TS_ASSERT_DELTA( phi_index, 4, TOLERATED_ERROR );
		TS_ASSERT_DELTA( theta_index, 33, TOLERATED_ERROR );

		Vector cc( -0.3425, -0.0166, 0.9394 ); // switch phi and theta
		core::scoring::sasa::get_legrand_orientation( cc, d, phi_index, theta_index, dist );
		TS_ASSERT_DELTA( phi_index, 4, TOLERATED_ERROR );
		TS_ASSERT_DELTA( theta_index, 33, TOLERATED_ERROR );

		//std::cout << "\n" << std::endl;


		// lines 5 of the database file sampling/SASA-angles.dat:
		// 6 6 6 7 7 7 7 7 7 7 7 5 5 5 5 5 18 18 18 18 18 18 18 17 17 17 17 17 17 28 28 28 28 28 28 28 27 27 27 27 27 27 38 38 38 38 38 38 38 37 37 37 37 37 45 45 45 45 45 45 45 45 6 6
		//
		// let's assume that the database is in theta,phi order
		// so line 5, or theta bin 5 (22.5-28.125 degrees, 0.392-0.4908 rad) and phi bin 63 (348.75-365.4degrees, 6.087-6.185radians) give a dot location of 6
		// use calc3d to get cartesian for ( 1.0, phi: 6.13, theta: 0.45 )
		dist = 1.0;
		Vector k( 0.4299, -0.0664, 0.90045);
		Vector l( 0, 0, 0 );
		core::scoring::sasa::get_legrand_orientation( k, l, phi_index, theta_index, dist );
		TS_ASSERT_DELTA( phi_index, 5, TOLERATED_ERROR );
		TS_ASSERT_DELTA( theta_index, 63, TOLERATED_ERROR );

		//std::cout << "\n" << std::endl;


		// I'm worried that the bottom 32 phi bins ( 180 to 360 degrees ) aren't being used correctly. That would be values
		// of phi and theta > pi (3.1415). Ok, let's pick a region of the database file to explore.
		// the last 3 lines
		// 28 27 27 27 27 27 27 27 27 27 27 27 27 37 37 37 37 37 37 37 37 37 37 37 37 37 6 6 6 6 6 6 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 17 17 17 17 17 17 17 17 17 17 17 17
		// 17 27 27 27 27 27 27 27 27 27 27 27 27 37 37 37 37 37 37 37 37 37 37 37 37 37 6 6 6 6 6 6 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 17 17 17 17 17 17 17 17 17 17 17 17
		// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

		// phi bin 63 (348.75-365.4degrees, 6.087-6.185radians), and theta bin 1 (0-5.625degrees, 0-0.098radians) is dot #17
		// using the calc3D program we can go from polar coordinates to cartesian. let's pick (1.0, phi: 6.14rad, theta: 0.05rad.

		dist = 1.00;
		Vector g( 0.04947, -0.007132, 0.9988 );
		Vector h( 0, 0, 0 );
		core::scoring::sasa::get_legrand_orientation( g, h, phi_index, theta_index, dist );
		TS_ASSERT_DELTA( phi_index, 63, TOLERATED_ERROR );
		TS_ASSERT_DELTA( theta_index, 1, TOLERATED_ERROR );

		//Vector gg( -0.1425, -0.007132, 0.9898 );
		//core::scoring::sasa::get_legrand_orientation( gg, h, phi_index, theta_index, dist ); // reverse phi/theta
		//TS_ASSERT_DELTA( phi_index, 63, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( theta_index, 1, TOLERATED_ERROR );

		//std::cout << "\n" << std::endl;


		// lines 60 of the database file sampling/SASA-angles.dat:
		// 28 28 28 28 26 26 26 26 26 26 38 38 38 38 38 38 38 36 36 36 36 36 36 45 45 45 45 45 45 45 4 4 4 4 4 7 7 7 7 7 7 7 3 3 3 3 3 3 18 18 18 18 18 18 18 16 16 16 16 16 16 28 28 28
		//
		// let's assume
		// so phi bin 60 (331.8-337.5 degrees, 5.791-5.890 rad) and theta bin 6 (28.125-33.75 degrees, 0.491-0.589radians) give a dot location of 26
		// use calc3d to get cartesian for ( 1.0, phi: 5.87, theta: 0.55 )
		dist = 1.0;
		Vector e( 0.4787, -0.2099, 0.8525 );
		Vector f( 0, 0, 0 );
		core::scoring::sasa::get_legrand( e, f, phi_index, theta_index, dist );
		TS_ASSERT_DELTA( phi_index, 60, TOLERATED_ERROR );
		TS_ASSERT_DELTA( theta_index, 6, TOLERATED_ERROR );

		//Vector ee( -0.3423, -0.2099, 0.9158 );
		//core::scoring::sasa::get_legrand_orientation( ee, f, phi_index, theta_index, dist );
		//TS_ASSERT_DELTA( phi_index, 60, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( theta_index, 6, TOLERATED_ERROR );

		//std::cout << "\n" << std::endl;


		// lines 60 of the database file sampling/SASA-angles.dat:
		// 28 28 28 28 26 26 26 26 26 26 38 38 38 38 38 38 38 36 36 36 36 36 36 45 45 45 45 45 45 45 4 4 4 4 4 7 7 7 7 7 7 7 3 3 3 3 3 3 18 18 18 18 18 18 18 16 16 16 16 16 16 28 28 28
		//
		// let's assume
		// so phi bin 60 (331.8-337.5 degrees, 5.791-5.890 rad) and theta bin 61 (337.5-343.125 degrees, 5.89-5.988 radians) give a dot location of 16
		// use calc3d to get cartesian for ( 1.0, phi: 5.87, theta: 5.93 ) same as ( 1.0, -0.41, -0.35 )
		dist = 1.0;
		Vector m( -0.3168, 0.1389, 0.9383 );
		Vector n( 0, 0, 0 );
		core::scoring::sasa::get_legrand_orientation( m, n, phi_index, theta_index, dist );
		TS_ASSERT_DELTA( phi_index, 60, TOLERATED_ERROR );
		TS_ASSERT_DELTA( theta_index, 61, TOLERATED_ERROR );

		//Vector mm( -0.3767, 0.1389, 0.9158 );
		//core::scoring::sasa::get_legrand_orientation( mm, n, phi_index, theta_index, dist );
		//TS_ASSERT_DELTA( phi_index, 60, TOLERATED_ERROR );
		//TS_ASSERT_DELTA( theta_index, 61, TOLERATED_ERROR );

		//std::cout << "\n" << std::endl;
	}


	void test_calc_total_sasa_glycine() {

		//chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		//core::pose::make_pose_from_sequence( pose, "G", *rsd_set );
		core::import_pose::pose_from_pdb( pose, "core/scoring/nonideal_glycine.pdb" );

		for ( Size ii=1; ii <= pose.n_residue(); ii+=3 ) {
			pose.set_phi( ii, -150.0 );
			pose.set_psi( ii, 150.0 );
			pose.set_omega( ii, 180.0 );
		}

		//ronj calling init_heavy_only leads to uninit'd values in the map
		//ronj what happens is that a vector1 exists for every position in the pose, and when init_heavy_only is called,
		//ronj the first heavy atom number of atoms at the beginning of the vector are set to 'true' or whatever the value
		//ronj type is and then the non-heavy atom indices of the vector are init'd with random values. This behavior occurs
		//ronj even though the vector at each residue position is resized to either the number of atoms or the number of heavy
		//ronj atoms depending on what function was called. It seems like the resize method is not really resizing the vector1
		//ronj at each position, leading to incorrect values being calculated in this method if it's used. Not about to try to
		//ronj debug a vector1 problem; therefore, I'm init'ing the values of the atom_subset atomid_map myself.
		//id::initialize_heavy_only( atom_subset, pose, true ); // this call leads to uninit'd values in the map, causing random results
		atom_subset.clear();
		atom_subset.resize( pose.n_residue() );
		for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
			atom_subset.resize( ii, pose.residue_type(ii).natoms(), false );
			for ( Size jj = 1; jj <= pose.residue_type(ii).nheavyatoms(); ++jj ) {
				atom_subset[ ii ][ jj ] = true;
			}
		}

		total_sasa = core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */, atom_subset );
		TS_ASSERT_DELTA( total_sasa, 223.8206, TOLERATED_ERROR );
		TS_ASSERT_DELTA( rsd_sasa[1], 223.8206, TOLERATED_ERROR );

		TS_ASSERT_DELTA( atom_sasa[ core::id::AtomID( 1, 1 )], 51.9289, TOLERATED_ERROR );
		TS_ASSERT_DELTA( atom_sasa[ core::id::AtomID( 2, 1 ) ], 86.9195, TOLERATED_ERROR );
		TS_ASSERT_DELTA( atom_sasa[ core::id::AtomID( 3, 1 ) ], 14.2736, TOLERATED_ERROR );

	}

	void test_calc_total_sasa_smallprotein() {

		pose = create_1ten_pdb_pose();
		core::pose::initialize_atomid_map( atom_subset, pose, true );

		total_sasa = core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */, atom_subset );
		TS_ASSERT_DELTA( total_sasa, 5240.4203, TOLERATED_ERROR );

		TS_ASSERT_DELTA( rsd_sasa[1], 66.9017, TOLERATED_ERROR );
		TS_ASSERT_DELTA( rsd_sasa[11], 117.2182, TOLERATED_ERROR );
		TS_ASSERT_DELTA( rsd_sasa[33], 0.0, TOLERATED_ERROR );

	}

	void test_calc_per_res_hydrophobic_sasa() {

		pose = create_1ten_pdb_pose();

		core::Real total_hydrophobic_sasa = 0.0;
		utility::vector1< core::Real > residue_sasa( pose.total_residue(), 0.0 );
		utility::vector1< core::Real > residue_hsasa( pose.total_residue(), 0.0 ); // hydrophobic SASA only
		total_hydrophobic_sasa = core::scoring::calc_per_res_hydrophobic_sasa( pose, residue_sasa, residue_hsasa, 1.4 /* probe radius */ );

		TS_ASSERT_DELTA( total_hydrophobic_sasa, 3212.4579, TOLERATED_ERROR );

		TS_ASSERT_DELTA( residue_hsasa[1], 27.7418, TOLERATED_ERROR );
		TS_ASSERT_DELTA( residue_hsasa[11], 62.4142, TOLERATED_ERROR );
		TS_ASSERT_DELTA( residue_hsasa[33], 0.0, TOLERATED_ERROR );

	}

};


