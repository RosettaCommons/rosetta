// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/pmut_scan/point_mut_scan.cxxtest.hh
/// @brief  test suite for the pilot app point_mut_scan
/// @author Ron Jacak (ron.jacak@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <basic/Tracer.hh>

// Unit headers
#include <protocols/pmut_scan/PointMutScanDriver.hh>
#include <protocols/pmut_scan/Mutant.hh>

#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

// Package Headers
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("test.protocols.pmut_scan.PointMutScanDriverTests");

// --------------- Test Class --------------- //

using namespace protocols::pmut_scan;

class PointMutScanDriverTests : public CxxTest::TestSuite {

public:

	bool suite_initialized;
	core::Real TOLERATED_ERROR;
	utility::vector1< core::pose::Pose > input_poses;

	utility::vector1< std::string > pdb_file_names;
	bool double_mutant_scan;
	bool output_mutant_structures;
	std::string list_file;


	// --------------- Suite-level Fixture --------------- //

	PointMutScanDriverTests() {
		suite_initialized = false;
	}

	virtual ~PointMutScanDriverTests() {}

	static PointMutScanDriverTests *createSuite() {
		return new PointMutScanDriverTests();
	}

	static void destroySuite( PointMutScanDriverTests *suite ) {
		delete suite;
	}

	void initialize_suite() {
		if ( suite_initialized ) return;
		suite_initialized = true;

		core_init_with_additional_options( "-mute core.pack.annealer core.pack.pack_rotamers core.pack.task core.pack.interaction_graph core.scoring core.io" );
		TOLERATED_ERROR = 0.1;

		// read in poses. only needed for 3 of the tests, but this is better than reading the PDBs in 3 times.
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/pmut_scan/shortloop.01.pdb.gz" , core::import_pose::PDB_file);
		input_poses.push_back( pose );

		core::import_pose::pose_from_file( pose, "protocols/pmut_scan/shortloop.02.pdb.gz" , core::import_pose::PDB_file);
		input_poses.push_back( pose );

		pdb_file_names.push_back( "protocols/pmut_scan/shortloop.01.pdb.gz" );
		pdb_file_names.push_back( "protocols/pmut_scan/shortloop.02.pdb.gz" );
	}


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case.

	// Shared initialization goes here.
	void setUp() {
		initialize_suite();

		double_mutant_scan = false;
		output_mutant_structures = false;
		list_file = "";
	}


	// Shared finalization goes here.
	void tearDown() {}


	// --------------- Test Cases --------------- //

	/// @brief tests the function fill_mutations_list for whether it reads in mutant list file properly. in particular,
	/// make sure single, double and higher order mutants are read in correctly.
	void test_fill_mutations_list() {
		TR << "Running test_fill_mutations_list..." << std::endl;

		list_file = "protocols/pmut_scan/mutants.txt";
		protocols::pmut_scan::PointMutScanDriver custom_test_driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );

		custom_test_driver.fill_mutations_list(); // reads in the file passed via the constructor

		TS_ASSERT_EQUALS( custom_test_driver.n_mutants(), core::Size(5) );

		// mutant #5
		//L S 30B T L N 30E W L T 31 Y

		MutationData md1( 'S', 'T', 2, 30, 'B', 'L' );
		MutationData md2( 'N', 'W', 5, 30, 'E', 'L' );
		MutationData md3( 'T', 'Y', 6, 31, ' ', 'L' );
		Mutant m;
		m.add_mutation( md1 );
		m.add_mutation( md2 );
		m.add_mutation( md3 );

		utility::vector1< protocols::pmut_scan::Mutant >::const_iterator it = custom_test_driver.mutants_begin();
		it += 4;
		TS_ASSERT( *it == m ); // also tests the operator== function!
	}

	/// @brief tests the function scan_for_mutations for whether single mutant scanning is set up correctly
	void test_fill_mutations_list_single_mutant_scan() {
		TR << "Running test_fill_mutations_list_single_mutant_scan..." << std::endl;

		protocols::pmut_scan::PointMutScanDriver driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );
		driver.fill_mutations_list();

		TS_ASSERT_EQUALS( driver.n_mutants(), core::Size(133) );

	}

	/// @brief tests the function scan_for_mutations for whether double mutant scanning is set up correctly
	void test_fill_mutations_list_double_mutant_scan() {
		TR << "Running test_fill_mutations_list_double_mutant_scan..." << std::endl;

		double_mutant_scan = true;
		protocols::pmut_scan::PointMutScanDriver custom_test_driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );

		custom_test_driver.fill_mutations_list();

		//1805 is the number of double mutants passing the distance filter (of 7581 total, which is (7 choose 2) * 19^2).  133 is the number of single mutants also included, now that I have changed double-mutant scan to include single mutants SML 11/2/12
		TS_ASSERT_EQUALS( custom_test_driver.n_mutants(), core::Size(1805 + 133) );

	}

	/// @brief tests the function are_residues_neighbors
	void test_calculate_neighbor_table() {
		TR << "Running test_calculate_neighbor_table..." << std::endl;

		protocols::pmut_scan::PointMutScanDriver driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );

		//core::graph::Graph neighbor_graph( input_poses[1].total_residue() );
		utility::vector1< utility::vector1< bool > > neighbors;
		driver.calculate_neighbor_table( input_poses[1], neighbors );

		TS_ASSERT( neighbors[ 1 ][ 2 ] );
		TS_ASSERT( !(neighbors[ 2 ][ 6 ]) );

	}

	/// @brief tests the function make_specific_mutant
	void test_make_specific_mutant() {
		TR << "Running test_make_specific_mutant..." << std::endl;

		protocols::pmut_scan::PointMutScanDriver driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );

		Mutant m;
		MutationData md1( 'T', 'Y', 6, 31, 'A', 'L' );
		m.add_mutation( md1 );

		// init some mutant poses
		utility::vector1< core::pose::Pose > mutant_poses( input_poses.size() );
		utility::vector1< core::pose::Pose > native_poses( input_poses.size() );
		for ( core::Size ii=1; ii <= input_poses.size(); ++ii ) {
			mutant_poses[ ii ] = input_poses[ ii ];
			native_poses[ ii ] = input_poses[ ii ];
		}

		driver.make_specific_mutant( mutant_poses, native_poses, m );

		TS_ASSERT_EQUALS( native_poses[1].residue( 1 ).name3(), "HIS" );
		TS_ASSERT_EQUALS( native_poses[1].residue( 6 ).name3(), "THR" );

		TS_ASSERT_EQUALS( mutant_poses[1].residue( 1 ).name3(), "HIS" );
		TS_ASSERT_EQUALS( mutant_poses[1].residue( 6 ).name3(), "TYR" );

	}

	/// @brief tests that the function make_specific_mutant for making mutants listed in a mutants.txt file is working
	void test_make_specific_mutant_insertion_code() {
		TR << "Running test_make_specific_mutant_insertion_code..." << std::endl;

		protocols::pmut_scan::PointMutScanDriver driver( pdb_file_names, double_mutant_scan, list_file, output_mutant_structures );

		MutationData md1( 'H', 'W', 1, 30, 'A', 'L' );
		MutationData md2( 'S', 'T', 2, 30, 'B', 'L' );
		MutationData md3( 'N', 'Q', 3, 30, 'C', 'L' );

		Mutant m;
		m.add_mutation( md1 );
		m.add_mutation( md2 );
		m.add_mutation( md3 );

		// init some mutant poses
		utility::vector1< core::pose::Pose > mutant_poses( input_poses.size() );
		utility::vector1< core::pose::Pose > native_poses( input_poses.size() );
		for ( core::Size ii=1; ii <= input_poses.size(); ++ii ) {
			mutant_poses[ ii ] = input_poses[ ii ];
			native_poses[ ii ] = input_poses[ ii ];
		}

		utility::vector1< Mutant >::const_iterator it;
		driver.make_specific_mutant( mutant_poses, native_poses, m );

		// assert that the native poses are not changed
		TS_ASSERT_EQUALS( native_poses[1].residue( 1 ).name3(), "HIS" );
		TS_ASSERT_EQUALS( native_poses[1].residue( 2 ).name3(), "SER" );
		TS_ASSERT_EQUALS( native_poses[1].residue( 3 ).name3(), "ASN" );

		TS_ASSERT_EQUALS( mutant_poses[1].residue( 1 ).name3(), "TRP" );
		TS_ASSERT_EQUALS( mutant_poses[1].residue( 2 ).name3(), "THR" );
		TS_ASSERT_EQUALS( mutant_poses[1].residue( 3 ).name3(), "GLN" );

	}


};
