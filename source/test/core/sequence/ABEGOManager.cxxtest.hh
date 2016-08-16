// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/ABEGOManager.cxxtest.hh
/// @brief  test suite for basic::ABEGOManager.cc
/// @author Colin A. Smith


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <core/sequence/ABEGOManager.hh>

#include <numeric/random/random.hh>
#include <utility/vector1.hh>
#include <string>

static basic::Tracer TR("test.core.sequence.ABEGOManager");

// --------------- Test Class --------------- //

class ABEGOManagerTest : public CxxTest::TestSuite {

public:

	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::sequence::ABEGOManager ABEGOManager;
	typedef core::pose::Pose Pose;

public:


	ABEGOManager abm;


public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
		// Being a smart pointer, g should be destructed and "free'd" correctly, but a g->delete_everything()
		// could be placed here, if desired.
	}

	// --------------- Test Cases --------------- //

	void test_abego() {

		//  ABEGO A( 'A', -180.0,   0.0,  -75.0,  50.0, false );
		//  ABEGO B( 'B', -180.0,   0.0,   50.0, 285.5, false );
		//  ABEGO E( 'E',    0.0, 180.5,  100.0, 260.5, false );
		//  ABEGO G( 'G',    0.0, 180.5, -100.0, 100.0, false );

		for ( Size ii=1; ii<=abm.total_number_abego(); ii++ ) {
			char symbol = abm.index2symbol( ii );
			TS_ASSERT( ii == abm.symbol2index( symbol ) );
		}

		// check omega cis
		char symbol = abm.index2symbol( abm.torsion2index_level2( 100.0, 100.0, 0 ) );
		TS_ASSERT( symbol == 'O' );
		TS_ASSERT( abm.check_rama( symbol, 100.0, 100.0, 0 ) );

		// check level_1
		for ( Size ii=0; ii<=72; ii++ ) {
			for ( Size jj=0; jj<=72; jj++ ) {
				Real phi = Real( ii*5 - 180.0 );
				Real psi = Real( jj*5 - 180.0 );
				Real omega = 180.0;
				char symbol = abm.index2symbol( abm.torsion2index_level1( phi, psi, omega ) );
				TS_ASSERT( abm.check_rama( symbol, phi, psi, omega ) );
			}
		}

		// check level_2
		for ( Size ii=0; ii<=72; ii++ ) {
			for ( Size jj=0; jj<=72; jj++ ) {
				Real phi = Real( ii*5 - 180.0 );
				Real psi = Real( jj*5 - 180.0 );
				Real omega = 180.0;
				char symbol = abm.index2symbol( abm.torsion2index_level2( phi, psi, omega ) );
				TS_ASSERT( abm.check_rama( symbol, phi, psi, omega ) );
			}
		}

		// check level_3
		for ( Size ii=0; ii<=72; ii++ ) {
			for ( Size jj=0; jj<=72; jj++ ) {
				Real phi = Real( ii*5 - 180.0 );
				Real psi = Real( jj*5 - 180.0 );
				Real omega = 180.0;
				char symbol = abm.index2symbol( abm.torsion2index_level3( phi, psi, omega ) );
				TS_ASSERT( abm.check_rama( symbol, phi, psi, omega ) );
			}
		}

		// check level_4
		for ( Size ii=0; ii<=72; ii++ ) {
			for ( Size jj=0; jj<=72; jj++ ) {
				Real phi = Real( ii*5 - 180.0 );
				Real psi = Real( jj*5 - 180.0 );
				Real omega = 180.0;
				char symbol = abm.index2symbol( abm.torsion2index_level4( phi, psi, omega ) );
				TS_ASSERT( abm.check_rama( symbol, phi, psi, omega ) );
			}
		}

		// check get_abego_string
		utility::vector1< String > abego;
		abego.push_back( "A" );
		abego.push_back( "B" );
		abego.push_back( "SP" );
		abego.push_back( "MN" );
		abego.push_back( "E" );
		abego.push_back( "O" );
		TS_ASSERT( "AB[SP][MN]EO" == abm.get_abego_string( abego ) );

		Pose pose;
		core::import_pose::pose_from_file( pose, "core/sequence/abego_test.pdb" , core::import_pose::PDB_file);
		// check get abego from pose
		TS_ASSERT( abm.get_abego_string( core::sequence::get_abego( pose, /* level */ 1 ) ) == "EBAAAGBO" );
		TS_ASSERT( abm.get_abego_string( core::sequence::get_abego( pose, /* level */ 3 ) ) == "EZAAAGSO" );
		TS_ASSERT( abm.get_abego_string( core::sequence::get_abego( pose, /* level */ 4 ) ) == "EZNNMGSO" );

	}

};


