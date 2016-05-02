// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/RotamerDots.cxxtest.hh
/// @brief  test suite for the rotamerdots and dotssphere classes
/// @author Ron Jacak

// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/id/AtomID.hh>
#include <core/pack/interaction_graph/RotamerDots.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/sasa.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

#include <core/import_pose/import_pose.hh>

// Utility Headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Test headers
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/ubyte.hh>


static basic::Tracer TR("test.core.pack.interaction_graph.rotamerdots");

using namespace core;
using namespace core::pack;
using namespace core::pack::interaction_graph;


// --------------- Test Class --------------- //

class RotamerDotsTests : public CxxTest::TestSuite {

public:

	// Shared data elements go here.
	pose::Pose pose;
	Size num_bytes;
	bool exclude_hydrogens;
	bool use_expanded_polar_atom_radii;

	// --------------- Suite-level Fixture --------------- //

	RotamerDotsTests() {

		// if the tests are run manually (or one suite at a time), that doesn't mute all of the tracer output by default.  Place
		// a mute here because the interaction graphs generate tons of debugging output (in DEBUG mode anyway).
		core_init_with_additional_options( "-no_optH -mute core.io core.mm" );

		// --- Pose ---
		// since this is a test suite, we don't want to read in PDB files from the command line.  just hardcode the tests to use
		// a predefined test PDB file
		core::import_pose::pose_from_file( pose, "core/pack/1l2y_renameH.pdb" , core::import_pose::PDB_file);

		num_bytes = 21;
		exclude_hydrogens = true;
		use_expanded_polar_atom_radii = true;

	}

	virtual ~RotamerDotsTests() {}

	static RotamerDotsTests *createSuite() {
		return new RotamerDotsTests();
	}

	static void destroySuite( RotamerDotsTests *suite ) {
		delete suite;
	}


	// --------------- Test Fixture --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case.

	void setUp() {}

	// Shared finalization goes here.
	// All memory allocated via OPs; objects should destroy themselves so nothing else to do here.
	void tearDown() {}


public:

	// --------------- Helper functions --------------- //

	void print_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) {
		for ( Size bb = 1; bb <= num_bytes; ++bb ) {
			if ( bb % 2 == 1 ) TR << (bb-1) * 8 << ":";
			for ( int index=7; index >= 0; index-- ) {
				TR << ( ( (int)values[ bb ] >> index ) & 1 );
			}
			TR << " ";
		}
		TR << std::endl;
	}

	// --------------- Test Cases --------------- //

	/// @details
	/// Tests to make sure when doing a design on only some residues that certain residues are indeed being treated and set
	/// as background nodes. If this array returns the wrong indices, background nodes are not being set properly.
	///
	void test_dotsphere_construction() {
		TR << "Running test_dotsphere_construction..." << std::endl;
		DotSphere ds;
		TS_ASSERT_EQUALS( ds.get_num_covered(), 0 ); // default constructed object should have zero for num_covered_

		TS_ASSERT_EQUALS( ds.get_dot_covered( 1 ), 0 ); // remember, dot ids are 1-based!
		TS_ASSERT_EQUALS( ds.get_dot_covered( 25 ), 0 );
		TS_ASSERT_EQUALS( ds.get_dot_covered( 162 ), 0 );
	}

	// @details
	// tests to make sure the member function get_total_dots is correct
	void test_dotsphere_get_total_dots() {
		TR << "Running test_dotsphere_get_total_dots..." << std::endl;
		DotSphere ds;
		TS_ASSERT_EQUALS( ds.get_total_dots(), (Size)162 );
	}

	void test_dotsphere_increment_count() {
		TR << "Running test_dotsphere_increment_count..." << std::endl;
		DotSphere ds;
		TS_ASSERT_EQUALS( ds.get_num_covered(), 0 ); // default constructed object should have zero for num_covered_

		utility::vector1< ObjexxFCL::ubyte > test_case_mask( num_bytes, ObjexxFCL::ubyte( 0 ) ); // create an empty mask vector

		// make up a fake masknum like 2549, which in the sampling/SASA-masks.dat file is
		// 97 128 119 254 255 63 6 0 136 255 255 39 0 192 127 7 0 0 105 0 0
		int masknum = 2549;
		for ( Size bb = 1; bb <= num_bytes; ++bb ) {
			test_case_mask[ bb ] |= core::scoring::get_masks()( bb, masknum );
		}
		//print_bit_string( test_case_mask );

		ds.increment_count( test_case_mask );
		//ds.print( std::cout );

		TS_ASSERT_EQUALS( ds.get_num_covered(), 71 ); // default constructed object should have zero for num_covered_

	}

	void test_dotsphere_get_num_covered() {
		TR << "Running test_dotsphere_get_num_covered..." << std::endl;
		DotSphere ds;
		utility::vector1< ObjexxFCL::ubyte > test_case_mask( num_bytes, ObjexxFCL::ubyte( 0 ) ); // create an empty mask vector

		for ( Size bb = 1; bb <= num_bytes; ++bb ) {
			test_case_mask[ bb ] |= core::scoring::get_masks()( bb, 1 );
		}
		ds.increment_count( test_case_mask );
		TS_ASSERT_EQUALS( ds.get_num_covered(), 0 ); // line 1 of the masks database file has all zeros in it


		for ( Size bb = 1; bb <= num_bytes; ++bb ) {
			test_case_mask[ bb ] |= core::scoring::get_masks()( bb, 202 );
		}
		//print_bit_string( test_case_mask );

		ds.increment_count( test_case_mask );
		//ds.print( std::cout );
		TS_ASSERT_EQUALS( ds.get_num_covered(), 1 ); // line 202 of the masks database file has a single 4 on it, the first byte

		TS_ASSERT_EQUALS( ds.get_num_uncovered(), 161 );

		// test a line that has all ones on it
		for ( Size bb = 1; bb <= num_bytes; ++bb ) {
			test_case_mask[ bb ] |= core::scoring::get_masks()( bb, 16000 );
		}
		//print_bit_string( test_case_mask );

		ds.increment_count( test_case_mask );
		//ds.print( std::cout );
		TS_ASSERT_EQUALS( ds.get_num_covered(), 161 );
		TS_ASSERT_EQUALS( ds.get_num_uncovered(), 1 );

		// test a random line
		utility::vector1< ObjexxFCL::ubyte > test_case_mask2( num_bytes, ObjexxFCL::ubyte( 0 ) ); // create an empty mask vector
		for ( Size bb = 1; bb <= num_bytes; ++bb ) {
			test_case_mask2[ bb ] |= core::scoring::get_masks()( bb, 3456 );
		}
		//print_bit_string( test_case_mask2 );

		DotSphere ds1;
		ds1.increment_count( test_case_mask2 );
		//ds1.print( std::cout );
		TS_ASSERT_EQUALS( ds1.get_num_covered(), 90 );

	}

	void test_dotsphere_operators() {
		TR << "Running test_dotsphere_operators..." << std::endl;
		DotSphere ds;
		utility::vector1< ObjexxFCL::ubyte > test_case_mask( num_bytes, ObjexxFCL::ubyte( 0 ) );
		for ( Size bb = 1; bb <= num_bytes; ++bb ) {
			test_case_mask[ bb ] |= core::scoring::get_masks()( bb, 2549 );
		}
		ds.increment_count( test_case_mask );
		TS_ASSERT_EQUALS( ds.get_dot_covered( 17 ), true );

		DotSphere ds1;
		ds1.increment_count( test_case_mask );
		TS_ASSERT_EQUALS( ds1.get_dot_covered( 17 ), true );
		ds -= ds1;
		TS_ASSERT_EQUALS( ds.get_dot_covered( 17 ), false ); // ds should be uncovered now
		TS_ASSERT_EQUALS( ds1.get_dot_covered( 17 ), true ); // but ds1 should still be covered

		ds1 += ds1;
		TS_ASSERT_EQUALS( ds1.get_dot_covered( 145 ), true );
		TS_ASSERT_EQUALS( ds1.get_dot_covered( 137 ), false );
		//ds1.print( std::cout );

	}

	void test_dotsphere_write_to_compact_array() {
		TR << "Running test_dotsphere_write_to_compact_array..." << std::endl;
		DotSphere ds;
		utility::vector1< ObjexxFCL::ubyte > test_case_mask( num_bytes, ObjexxFCL::ubyte( 0 ) );
		utility::vector1< ObjexxFCL::ubyte > compact_array( num_bytes, ObjexxFCL::ubyte( 0 ) );

		for ( Size bb = 1; bb <= num_bytes; ++bb ) {
			test_case_mask[ bb ] |= core::scoring::get_masks()( bb, 202 );
		}
		for ( Size ii = 1; ii <= 10; ++ii ) {
			ds.increment_count( test_case_mask );
		}
		//ds.print( std::cout );
		ds.write_to_compact_array( compact_array );
		TS_ASSERT_EQUALS( compact_array[ 1 ], 0x04 );
	}

	void test_rotamerdots_construction() {
		TR << "Running test_rotamerdots_construction..." << std::endl;
		RotamerDots rd;
		TS_ASSERT_EQUALS( rd.state_unassigned(), true );

		conformation::ResidueOP res = pose.residue(1).clone();
		RotamerDots rd1( res );
		TS_ASSERT_EQUALS( rd1.state_unassigned(), false );
		TS_ASSERT_EQUALS( rd1.get_num_atoms(), 16 );
		TS_ASSERT_EQUALS( rd1.get_atom_radius(1), 1.50 );
		numeric::xyzVector< Real > xyz( -8.901, 4.127, -0.555 );
		TS_ASSERT_DELTA( rd1.get_atom_coords_xyz( 1 ).x(), xyz.x(), 0.01 );
		TS_ASSERT_DELTA( rd1.get_atom_coords_xyz( 1 ).y(), xyz.y(), 0.01 );
		TS_ASSERT_DELTA( rd1.get_atom_coords_xyz( 1 ).z(), xyz.z(), 0.01 );

		RotamerDots rd2 = rd1; // tests both operator= and copy() method
		TS_ASSERT_EQUALS( rd2.state_unassigned(), false );
		TS_ASSERT_EQUALS( rd2.get_num_atoms(), 16 );
		TS_ASSERT_EQUALS( rd2.get_atom_radius(1), 1.50 );
		TS_ASSERT_DELTA( rd2.get_atom_coords_xyz( 1 ).x(), xyz.x(), 0.01 );
		TS_ASSERT_DELTA( rd2.get_atom_coords_xyz( 1 ).y(), xyz.y(), 0.01 );
		TS_ASSERT_DELTA( rd2.get_atom_coords_xyz( 1 ).z(), xyz.z(), 0.01 );

		RotamerDots rd3(rd1); // tests copy constructor
		TS_ASSERT_EQUALS( rd3.state_unassigned(), false );
		TS_ASSERT_EQUALS( rd3.get_num_atoms(), 16 );
		TS_ASSERT_EQUALS( rd3.get_atom_radius(1), 1.50 );
		TS_ASSERT_DELTA( rd3.get_atom_coords_xyz( 1 ).x(), xyz.x(), 0.01 );

		// exclude hydrogens
		RotamerDots rd4( res, exclude_hydrogens ); // tests constructor with extra bool parameter
		TS_ASSERT_EQUALS( rd4.state_unassigned(), false );
		TS_ASSERT_EQUALS( rd4.get_num_atoms(), 8 );
		TS_ASSERT_EQUALS( rd4.get_atom_radius(1), 1.50 );

		RotamerDots rd5 = rd4; // tests operator= in case of excluded hydrogens
		TS_ASSERT_EQUALS( rd5.state_unassigned(), false );
		TS_ASSERT_EQUALS( rd5.get_num_atoms(), 8 );
		TS_ASSERT_EQUALS( rd5.get_atom_radius(1), 1.50 );

		// exclude hydrogens and expanded polars
		RotamerDots rd6( res, exclude_hydrogens, use_expanded_polar_atom_radii ); // tests constructor with 2 extra bool parameters
		TS_ASSERT_EQUALS( rd6.state_unassigned(), false );
		TS_ASSERT_EQUALS( rd6.get_num_atoms(), 8 );
		TS_ASSERT_EQUALS( rd6.get_atom_radius(1), 2.50 ); // first atom is Nlys, so this should return a different value from before

		RotamerDots rd7 = rd6; // tests operator= in case of excluded hydrogens and expanded polar sasa
		TS_ASSERT_EQUALS( rd7.state_unassigned(), false );
		TS_ASSERT_EQUALS( rd7.get_num_atoms(), 8 );
		TS_ASSERT_EQUALS( rd7.get_atom_radius(1), 2.50 );

	}


	void test_rotamerdots_overlaps() {
		TR << "Running test_rotamerdots_overlaps..." << std::endl;
		conformation::ResidueOP res1 = pose.residue(1).clone();
		conformation::ResidueOP res2 = pose.residue(2).clone();
		conformation::ResidueOP res10 = pose.residue(10).clone();
		conformation::ResidueOP res5 = pose.residue(5).clone();
		conformation::ResidueOP res16 = pose.residue(16).clone();

		RotamerDots rd1( res1 );
		RotamerDots rd2( res2 );
		TS_ASSERT_EQUALS( rd1.overlaps( rd2 ), true );

		RotamerDots rd10( res10 );
		TS_ASSERT_EQUALS( rd1.overlaps( rd10 ), false );

		// do all of the same above except exclude hydrogens this time
		RotamerDots rd3( res1, exclude_hydrogens );
		RotamerDots rd4( res2, exclude_hydrogens );
		TS_ASSERT_EQUALS( rd3.overlaps( rd4 ), true );

		RotamerDots rd10b( res10, exclude_hydrogens );
		TS_ASSERT_EQUALS( rd3.overlaps( rd10b ), false );

		// and then exclude hydrogens and expand polars
		RotamerDots rd5( res5, exclude_hydrogens, use_expanded_polar_atom_radii );
		RotamerDots rd16( res16, exclude_hydrogens, use_expanded_polar_atom_radii );
		TS_ASSERT_EQUALS( rd5.overlaps( rd16 ), true );
		TS_ASSERT_EQUALS( rd3.overlaps( rd16 ), false );

	}


	void test_rotamerdots_increment_self_overlap() {
		TR << "Running test_rotamerdots_increment_self_overlap..." << std::endl;
		conformation::ResidueOP res = pose.residue(1).clone();
		RotamerDots rd1( res );

		//rd1.print( std::cout );
		rd1.increment_self_overlap();
		//rd1.print( std::cout );

		DotSphere ds1 = rd1.get_atom_counts()[ 1 ];
		TS_ASSERT_EQUALS( ds1.get_num_covered(), 149 );

		DotSphere ds16 = rd1.get_atom_counts()[ 16 ];
		TS_ASSERT_EQUALS( ds16.get_num_uncovered(), 35 );
		TS_ASSERT_EQUALS( rd1.get_num_uncovered( 16 ), 35 );

		// do all of the same above except exclude hydrogens this time
		RotamerDots rd2( res, exclude_hydrogens );

		//rd2.print( std::cout );
		rd2.increment_self_overlap();
		//rd2.print( std::cout );

		ds1 = rd2.get_atom_counts()[ 1 ];
		TS_ASSERT_EQUALS( ds1.get_num_covered(), 102 );

		DotSphere ds8 = rd2.get_atom_counts()[ 8 ];
		TS_ASSERT_EQUALS( ds8.get_num_uncovered(), 83 );
		TS_ASSERT_EQUALS( rd2.get_num_uncovered( 8 ), 83 );

		// and finally exclude hydrogens and expand polars
		RotamerDots rd3( res, exclude_hydrogens, use_expanded_polar_atom_radii );
		//rd3.print( std::cout );
		rd3.increment_self_overlap();
		//rd3.print( std::cout );

		ds1 = rd3.get_atom_counts()[ 1 ];
		TS_ASSERT_EQUALS( ds1.get_num_covered(), 89 ); // this will now be counting only the expanded counts
		ds8 = rd3.get_atom_counts()[ 8 ];
		TS_ASSERT_EQUALS( ds8.get_num_uncovered(), 87 ); // this will now be counting only the expanded counts
		TS_ASSERT_EQUALS( rd3.get_num_uncovered( 8 ), 87 ); // this will now be counting only the expanded counts

	}

	void test_rotamerdots_increment_this_and_cache() {
		TR << "Running test_rotamerdots_increment_this_and_cache..." << std::endl;
		conformation::ResidueOP res1 = pose.residue(1).clone();
		conformation::ResidueOP res2 = pose.residue(2).clone();
		RotamerDots rd1( res1 );
		RotamerDots rd2( res2 );
		utility::vector1< utility::vector1< bool > > atom_atom_overlap( res2->natoms(), utility::vector1< bool >( res1->natoms(), false ) );

		RotamerDotsCache rdc2( res2->natoms() );
		rd1.increment_this_and_cache( rd2, rdc2, atom_atom_overlap );

		//rd1.print( std::cout );
		//rd2.print( std::cout );

		utility::vector1< DotSphere > const & v_ds = rd1.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 51 );
		TS_ASSERT_EQUALS( v_ds[ 5 ].get_num_covered(), 39 );
		TS_ASSERT_EQUALS( v_ds[ 16 ].get_num_covered(), 0 );

		utility::vector1< DotSphere > const & v_ds2 = rd2.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds2[ 1 ].get_num_covered(), 0 ); // res2 is not altered in increment_this_and_cache
		TS_ASSERT_EQUALS( v_ds2[ 19 ].get_num_covered(), 0 );

		// do all of the same above except exclude hydrogens this time
		RotamerDots rd3( res1, exclude_hydrogens );
		RotamerDots rd4( res2, exclude_hydrogens );
		utility::vector1< utility::vector1< bool > > atom_atom_overlap_3( res2->nheavyatoms(), utility::vector1< bool >( res1->nheavyatoms(), false ) );

		RotamerDotsCache rdc4( res2->nheavyatoms() );
		rd3.increment_this_and_cache( rd4, rdc4, atom_atom_overlap_3 );

		//rd3.print( std::cout );
		//rd4.print( std::cout );

		utility::vector1< DotSphere > const & v_ds3 = rd3.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds3[ 1 ].get_num_covered(), 45 );
		TS_ASSERT_EQUALS( v_ds3[ 5 ].get_num_covered(), 38 );
		TS_ASSERT_EQUALS( v_ds3[ 8 ].get_num_covered(), 7 );

		utility::vector1< DotSphere > const & v_ds4 = rd4.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds4[ 1 ].get_num_covered(), 0 ); // res2 is not altered in increment_this_and_cache
		TS_ASSERT_EQUALS( v_ds4[ 8 ].get_num_covered(), 0 );
	}

	void test_rotamerdots_increment_this_and_cache_expanded_polars() {
		TR << "Running test_rotamerdots_increment_this_and_cache_expanded_polars..." << std::endl;
		conformation::ResidueOP res1 = pose.residue(1).clone();
		conformation::ResidueOP res2 = pose.residue(2).clone();
		RotamerDots rd5( res1, exclude_hydrogens, use_expanded_polar_atom_radii );
		RotamerDots rd6( res2, exclude_hydrogens, use_expanded_polar_atom_radii );
		utility::vector1< utility::vector1< bool > > atom_atom_overlap_1( res2->natoms(), utility::vector1< bool >( res1->natoms(), false ) );

		//rd5.print( std::cout );

		RotamerDotsCache rdc5( res2->nheavyatoms() );
		rd5.increment_this_and_cache( rd6, rdc5, atom_atom_overlap_1 );

		//rd5.print( std::cout );
		//rd6.print( std::cout );

		utility::vector1< DotSphere > const & v_ds5 = rd5.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds5[ 1 ].get_num_covered(), 54 );
		TS_ASSERT_EQUALS( v_ds5[ 5 ].get_num_covered(), 66 );
		TS_ASSERT_EQUALS( v_ds5[ 8 ].get_num_covered(), 21 );

		utility::vector1< DotSphere > const & v_ds6 = rd6.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds6[ 1 ].get_num_covered(), 0 ); // res2 is not altered in increment_this_and_cache
		TS_ASSERT_EQUALS( v_ds6[ 8 ].get_num_covered(), 0 );

		/*TR << "test_rotamerdots_increment_this_and_cache_expanded_polars(): rd5 expanded polar atom sasas: [ ";
		for ( Size ii=1; ii <= res1->nheavyatoms(); ++ii ) {
		TR << rd5.get_atom_sasa( ii ) << ", ";
		}
		TR << "]" << std::endl;*/

		TS_ASSERT_DELTA( rd5.get_atom_sasa(1), 127.4229, 0.1 );
		TS_ASSERT_DELTA( rd5.get_atom_sasa(5), 79.6272, 0.1 );
		TS_ASSERT_DELTA( rd5.get_atom_sasa(8), 179.4006, 0.1 );
	}

	void test_rotamerdots_increment_and_decrement_from_cached() {
		TR << "Running test_rotamerdots_increment_and_decrement_from_cached..." << std::endl;
		conformation::ResidueOP res1 = pose.residue(1).clone();
		RotamerDots rd1( res1 );
		//rd1.print( std::cout );
		utility::vector1< DotSphere > const & v_ds = rd1.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 0 );

		RotamerDotsCache rdc1( res1->natoms() );
		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > test_case_mask( res1->natoms(), utility::vector1< ObjexxFCL::ubyte >( num_bytes, ObjexxFCL::ubyte( 0 ) ) );
		for ( Size aa = 1; aa <= test_case_mask.size(); ++aa ) {
			for ( Size bb = 1; bb <= num_bytes; ++bb ) {
				test_case_mask[ aa ][ bb ] |= core::scoring::get_masks()( bb, 1103 );
			}
		}
		rdc1.increment_count( test_case_mask );
		rd1.increment_from_cached( rdc1 );
		//rd1.print( std::cout );
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 2 );
		TS_ASSERT_EQUALS( v_ds[ 16 ].get_num_covered(), 2 );

		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > test_case_mask2( res1->natoms(), utility::vector1< ObjexxFCL::ubyte >( num_bytes, ObjexxFCL::ubyte( 0 ) ) );
		for ( Size aa = 1; aa <= test_case_mask2.size(); ++aa ) {
			for ( Size bb = 1; bb <= num_bytes; ++bb ) {
				test_case_mask2[ aa ][ bb ] |= core::scoring::get_masks()( bb, 1102 );
			}
		}
		RotamerDotsCache rdc2( res1->natoms() );
		rdc2.increment_count( test_case_mask2 );
		rd1.decrement_from_cached( rdc2 );
		//rd1.print( std::cout );
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 1 );
		TS_ASSERT_EQUALS( v_ds[ 16 ].get_num_covered(), 1 );
	}

	void test_rotamerdots_increment_and_decrement_from_cached_excl_hydrogens() {
		TR << "Running test_rotamerdots_increment_and_decrement_from_cached_excl_hydrogens..." << std::endl;
		// same as the test above except exclude hydrogens this time

		conformation::ResidueOP res1 = pose.residue(1).clone();
		RotamerDots rd1( res1, exclude_hydrogens );
		//rd1.print( std::cout );
		utility::vector1< DotSphere > const & v_ds = rd1.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 0 );

		RotamerDotsCache rdc1( res1->nheavyatoms() );
		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > test_case_mask( res1->nheavyatoms(), utility::vector1< ObjexxFCL::ubyte >( num_bytes, ObjexxFCL::ubyte( 0 ) ) );
		for ( Size aa = 1; aa <= test_case_mask.size(); ++aa ) {
			for ( Size bb = 1; bb <= num_bytes; ++bb ) {
				test_case_mask[ aa ][ bb ] |= core::scoring::get_masks()( bb, 1103 );
			}
		}
		rdc1.increment_count( test_case_mask );
		rd1.increment_from_cached( rdc1 );
		//rd1.print( std::cout );
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 2 );
		TS_ASSERT_EQUALS( v_ds[ 8 ].get_num_covered(), 2 );

		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > test_case_mask2( res1->nheavyatoms(), utility::vector1< ObjexxFCL::ubyte >( num_bytes, ObjexxFCL::ubyte( 0 ) ) );
		for ( Size aa = 1; aa <= test_case_mask2.size(); ++aa ) {
			for ( Size bb = 1; bb <= num_bytes; ++bb ) {
				test_case_mask2[ aa ][ bb ] |= core::scoring::get_masks()( bb, 1102 );
			}
		}
		RotamerDotsCache rdc2( res1->nheavyatoms() );
		rdc2.increment_count( test_case_mask2 );
		rd1.decrement_from_cached( rdc2 );
		//rd1.print( std::cout );
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 1 );
		TS_ASSERT_EQUALS( v_ds[ 8 ].get_num_covered(), 1 );

	}

	void test_rotamerdots_increment_and_decrement_from_cached_expanded_polars() {
		TR << "Running test_rotamerdots_increment_and_decrement_from_cached_expanded_polars..." << std::endl;
		conformation::ResidueOP res1 = pose.residue(1).clone();

		// use a more complicated overlap mask and also keep expanded polar atom sasas
		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > test_case_mask( res1->nheavyatoms(), utility::vector1< ObjexxFCL::ubyte >( num_bytes, ObjexxFCL::ubyte( 0 ) ) );
		for ( Size aa = 1; aa <= test_case_mask.size(); ++aa ) {
			for ( Size bb = 1; bb <= num_bytes; ++bb ) {
				test_case_mask[ aa ][ bb ] |= core::scoring::get_masks()( bb, 3456 );
			}
		}

		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > test_case_mask2( res1->nheavyatoms(), utility::vector1< ObjexxFCL::ubyte >( num_bytes, ObjexxFCL::ubyte( 0 ) ) );
		for ( Size aa = 1; aa <= test_case_mask2.size(); ++aa ) {
			for ( Size bb = 1; bb <= num_bytes; ++bb ) {
				test_case_mask2[ aa ][ bb ] |= core::scoring::get_masks()( bb, 1102 );
			}
		}

		RotamerDotsCache rdc1( res1->nheavyatoms() );
		rdc1.increment_count( test_case_mask ); // no _ep specific version of this method anymore!
		RotamerDots rd1( res1, exclude_hydrogens, use_expanded_polar_atom_radii );
		//rd1.print( std::cout );

		rd1.increment_from_cached( rdc1 );
		//rd1.print( std::cout );

		RotamerDotsCache rdc2( res1->nheavyatoms() );
		rdc2.increment_count( test_case_mask2 );
		rd1.decrement_from_cached( rdc2 );
		//rd1.print( std::cout );

		utility::vector1< DotSphere > const & v_ds = rd1.get_atom_counts(); // will return counts only for expanded polar SASA
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 89 );
		TS_ASSERT_EQUALS( v_ds[ 8 ].get_num_covered(), 89 );

		/*TR << "test_rotamerdots_increment_and_decrement_from_cached_expanded_polars(): rd1 expanded polar atom sasas: [ ";
		for ( Size ii=1; ii <= res1->nheavyatoms(); ++ii ) {
		TR << rd1.get_atom_sasa( ii ) << ", ";
		}
		TR << "]" << std::endl;*/

		TS_ASSERT_DELTA( rd1.get_atom_sasa(1), 86.1285, 0.1 );
		TS_ASSERT_DELTA( rd1.get_atom_sasa(8), 92.8811, 0.1 );
	}

	void test_rotamerdots_increment_both_and_cache() {
		TR << "Running test_rotamerdots_increment_both_and_cache..." << std::endl;
		conformation::ResidueOP res1 = pose.residue(1).clone();
		conformation::ResidueOP res2 = pose.residue(2).clone();
		RotamerDots rd1( res1 );
		RotamerDots rd2( res2 );
		utility::vector1< utility::vector1< bool > > atom_atom_overlap( res2->natoms(), utility::vector1< bool >( res1->natoms(), false ) );

		RotamerDotsCache rdc2( res2->natoms() );
		RotamerDotsCache rdc1( res1->natoms() );
		rd1.increment_both_and_cache( rd2, rdc2, rdc1, atom_atom_overlap );

		utility::vector1< DotSphere > const & v_ds = rd1.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 51 );
		TS_ASSERT_EQUALS( v_ds[ 5 ].get_num_covered(), 39 );
		TS_ASSERT_EQUALS( v_ds[ 16 ].get_num_covered(), 0 );

		utility::vector1< DotSphere > const & v_ds2 = rd2.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds2[ 1 ].get_num_covered(), 89 );
		TS_ASSERT_EQUALS( v_ds2[ 10 ].get_num_covered(), 57 );
		TS_ASSERT_EQUALS( v_ds2[ 19 ].get_num_covered(), 2 );

		//rd1.print( std::cout );
		//rd2.print( std::cout );
	}

	void test_rotamerdots_increment_both_and_cache_excl_hydrogens() {
		TR << "Running test_rotamerdots_increment_both_and_cache_excl_hydrogens..." << std::endl;
		conformation::ResidueOP res1 = pose.residue(1).clone();
		conformation::ResidueOP res2 = pose.residue(2).clone();
		RotamerDots rd1( res1, exclude_hydrogens );
		RotamerDots rd2( res2, exclude_hydrogens );
		utility::vector1< utility::vector1< bool > > atom_atom_overlap( res2->nheavyatoms(), utility::vector1< bool >( res1->nheavyatoms(), false ) );

		RotamerDotsCache rdc2( res2->nheavyatoms() );
		RotamerDotsCache rdc1( res1->nheavyatoms() );
		rd1.increment_both_and_cache( rd2, rdc2, rdc1, atom_atom_overlap );

		utility::vector1< DotSphere > const & v_ds1 = rd1.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds1[ 1 ].get_num_covered(), 45 );
		TS_ASSERT_EQUALS( v_ds1[ 5 ].get_num_covered(), 38 );
		TS_ASSERT_EQUALS( v_ds1[ 8 ].get_num_covered(), 7 );

		utility::vector1< DotSphere > const & v_ds2 = rd2.get_atom_counts();
		TS_ASSERT_EQUALS( v_ds2[ 1 ].get_num_covered(), 88 );
		TS_ASSERT_EQUALS( v_ds2[ 4 ].get_num_covered(), 36 );
		TS_ASSERT_EQUALS( v_ds2[ 7 ].get_num_covered(), 13 );

		//rd1.print( std::cout );
		//rd2.print( std::cout );
	}

	void test_rotamerdots_increment_both_and_cache_expanded_polars() {
		TR << "Running test_rotamerdots_increment_both_and_cache_expanded_polars..." << std::endl;
		conformation::ResidueOP res1 = pose.residue(1).clone();
		conformation::ResidueOP res2 = pose.residue(2).clone();
		RotamerDots rd1( res1, exclude_hydrogens, use_expanded_polar_atom_radii );
		RotamerDots rd2( res2, exclude_hydrogens, use_expanded_polar_atom_radii );
		utility::vector1< utility::vector1< bool > > atom_atom_overlap( res2->nheavyatoms(), utility::vector1< bool >( res1->nheavyatoms(), false ) );

		RotamerDotsCache rdc2( res2->nheavyatoms() );
		RotamerDotsCache rdc1( res1->nheavyatoms() );
		rd1.increment_both_and_cache( rd2, rdc2, rdc1, atom_atom_overlap );

		utility::vector1< DotSphere > const & v_ds1 = rd1.get_atom_counts(); // no _ep specific version of this method anymore!
		TS_ASSERT_EQUALS( v_ds1[ 1 ].get_num_covered(), 54 );
		TS_ASSERT_EQUALS( v_ds1[ 5 ].get_num_covered(), 66 );
		TS_ASSERT_EQUALS( v_ds1[ 8 ].get_num_covered(), 21 );

		utility::vector1< DotSphere > const & v_ds2 = rd2.get_atom_counts(); // no _ep specific version of this method anymore!
		TS_ASSERT_EQUALS( v_ds2[ 1 ].get_num_covered(), 85 );
		TS_ASSERT_EQUALS( v_ds2[ 4 ].get_num_covered(), 46 );
		TS_ASSERT_EQUALS( v_ds2[ 7 ].get_num_covered(), 30 );

		//rd1.print( std::cout );
		//rd2.print( std::cout );

		/*TR << "test_rotamerdots_increment_both_and_cache_expanded_polars(): rd1 (expanded polar) atom sasas: [ ";
		for ( Size ii=1; ii <= res1->nheavyatoms(); ++ii ) {
		TR << rd1.get_atom_sasa( ii ) << ", ";
		}
		TR << "]" << std::endl;
		TR << "test_rotamerdots_increment_both_and_cache_expanded_polars(): rd2 (expanded polar) atom sasas: [ ";
		for ( Size ii=1; ii <= res2->nheavyatoms(); ++ii ) {
		TR << rd2.get_atom_sasa( ii ) << ", ";
		}
		TR << "]" << std::endl;*/

		TS_ASSERT_DELTA( rd1.get_atom_sasa(1), 127.4229, 0.1 );
		TS_ASSERT_DELTA( rd1.get_atom_sasa(5), 79.6272, 0.1 );
		TS_ASSERT_DELTA( rd1.get_atom_sasa(8), 179.4006, 0.1 );

		TS_ASSERT_DELTA( rd2.get_atom_sasa(1), 97.9705, 0.1 );
		TS_ASSERT_DELTA( rd2.get_atom_sasa(4), 129.9331, 0.1 );
		TS_ASSERT_DELTA( rd2.get_atom_sasa(7), 109.4874, 0.1 );
	}

	void test_rotamerdotscache_construction() {
		TR << "Running test_rotamerdotscache_construction..." << std::endl;
		conformation::ResidueOP res = pose.residue(1).clone();

		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > test_case_mask( res->natoms(), utility::vector1< ObjexxFCL::ubyte >( num_bytes, ObjexxFCL::ubyte( 0 ) ) );
		// store that line in the database file 16 times over in a vector
		for ( Size aa = 1; aa <= test_case_mask.size(); ++aa ) {
			for ( Size bb = 1; bb <= num_bytes; ++bb ) {
				test_case_mask[ aa ][ bb ] |= core::scoring::get_masks()( bb, 2549 );
			}
		}

		RotamerDotsCache rdc1( res->natoms() );
		rdc1.increment_count( test_case_mask );
		utility::vector1< DotSphere > const & v_ds = rdc1.get_atom_counts();
		//for ( Size ii=1; ii <= v_ds.size(); ++ii ) { v_ds[ ii ].print( std::cout ); }

		TS_ASSERT_EQUALS( v_ds[ 1 ].get_num_covered(), 71 );

		// do all of the same above except exclude hydrogens this time
		RotamerDotsCache rdc2( res->nheavyatoms() );
		test_case_mask.resize( res->nheavyatoms() ); // will truncate some data from the vector
		rdc2.increment_count( test_case_mask );
		utility::vector1< DotSphere > const & v_ds2 = rdc2.get_atom_counts();
		//for ( Size ii=1; ii <= v_ds2.size(); ++ii ) { v_ds2[ ii ].print( std::cout ); }

		TS_ASSERT_EQUALS( v_ds2[ 1 ].get_num_covered(), 71 );
		TS_ASSERT_EQUALS( v_ds2[ 8 ].get_num_covered(), 71 );


		// and finally keep polar atoms - though, there's not really a good way to test that polar sasa is being kept correctly
		// using the RDC API.
		RotamerDotsCache rdc3( res->nheavyatoms() );
		rdc3.increment_count( test_case_mask ); // no _ep specific version of this method anymore!
		utility::vector1< DotSphere > const & v_ds3 = rdc3.get_atom_counts();
		//for ( Size ii=1; ii <= v_ds3.size(); ++ii ) { v_ds3[ ii ].print( std::cout ); }

		TS_ASSERT_EQUALS( v_ds3[ 1 ].get_num_covered(), 71 );
		TS_ASSERT_EQUALS( v_ds3[ 5 ].get_num_covered(), 71 );
		TS_ASSERT_EQUALS( v_ds3[ 8 ].get_num_covered(), 71 );
	}

	void test_invrotamerdots() {
		TR << "Running test_invrotamerdots..." << std::endl;
		using namespace core::id;

		RotamerDots rd1( pose.residue(1).clone(), true, true );
		RotamerDots rd2( pose.residue(2).clone(), true, true );
		RotamerDots rd3( pose.residue(3).clone(), true, true );
		RotamerDots rd4( pose.residue(2).clone(), true, true );
		rd1.increment_self_overlap();
		rd2.increment_self_overlap();
		rd1.increment_both( rd2 );
		InvRotamerDots inv1, inv2;
		inv1.setup_from_rotamer_dots( rd1 );
		inv2.setup_from_rotamer_dots( rd2 );

		for ( Size ii = 1; ii <= pose.residue_type(1).nheavyatoms(); ++ii ) {
			for ( Size jj = 1; jj <= 162; ++jj ) { /// ugh... fix this constant, andrew.
				TS_ASSERT_EQUALS( rd1.get_atom_counts()[ ii ].get_dot_covered( jj ), ! inv1.dot_exposed( ii, jj ) );
				if ( rd1.get_atom_counts()[ ii ].get_dot_covered( jj ) == inv1.dot_exposed( ii, jj ) ) {
					std::cout << "Error on residue 1 on atom " << ii << " " << jj << " " << rd1.get_atom_counts()[ ii ].get_dot_covered( jj ) << " vs " << inv1.dot_exposed( ii, jj )  << std::endl;
				}
			}
		}

		for ( Size ii = 1; ii <= pose.residue_type(2).nheavyatoms(); ++ii ) {
			for ( Size jj = 1; jj <= 162; ++jj ) { /// ugh... fix this constant, andrew.
				TS_ASSERT_EQUALS( rd2.get_atom_counts()[ ii ].get_dot_covered( jj ), ! inv2.dot_exposed( ii, jj ) );
				if ( rd2.get_atom_counts()[ ii ].get_dot_covered( jj ) == inv2.dot_exposed( ii, jj ) ) {
					std::cout << "Error on residue 2 on atom " << ii << " " << jj << " " << rd2.get_atom_counts()[ ii ].get_dot_covered( jj ) << " vs " << inv2.dot_exposed( ii, jj )  << std::endl;
				}
			}
		}
	}

	void test_invrotamerdots_overlap_is_buried() {
		TR << "Running test_invrotamerdots_overlap_is_buried..." << std::endl;
		using namespace core::id;

		utility::vector1< RotamerDotsOP > rdots( 20 );
		utility::vector1< InvRotamerDotsOP > invdots( 20 );
		for ( Size ii = 1; ii <= 20; ++ii ) {
			rdots[ ii ] = RotamerDotsOP( new RotamerDots( pose.residue(ii).clone(), true, true ) );
			rdots[ ii ]->increment_self_overlap();
		}
		for ( Size ii = 1; ii <= 20; ++ii ) {
			for ( Size jj = ii+1; jj <= 20; ++jj ) {
				rdots[ii]->increment_both( *rdots[jj] );
			}
			invdots[ii] = InvRotamerDotsOP( new InvRotamerDots() );
			invdots[ii]->setup_from_rotamer_dots( *rdots[ii] );
		}
		//invdots[17]->write_exposed_dots_to_kinemage( std::cout, true );

		AtomID asn1_cb( pose.residue_type(1).atom_index( "CB" ), 1 );
		AtomID leu2_cb( pose.residue_type(2).atom_index( "CB" ), 2 );
		AtomID leu2_cd2( pose.residue_type(2).atom_index( "CD2" ), 2 );
		AtomID ile4_cg2( pose.residue_type(4).atom_index( "CG2" ), 4 );
		AtomID lys8_cb( pose.residue_type(8).atom_index( "CB" ), 8 );

		AtomID pro17_cb( pose.residue_type(17).atom_index( "CB" ), 17 );
		AtomID pro17_cd( pose.residue_type(17).atom_index( "CD" ), 17 );
		AtomID pro17_ca( pose.residue_type(17).atom_index( "CA" ), 17 );

		// test that exposed overlap is indeed called exposed
		TS_ASSERT( invdots[leu2_cb.rsd()]->atom_overlap_is_exposed( leu2_cb.atomno(), leu2_cd2.atomno() ));

		// test that buried overlap is called buried? although in PyMOL this looks very much buried the function is returning
		// true because both atoms have exposed dots, and it just so happens that there are dots adjacent to the intesection
		// circle on both atoms. so this is really testing that exposed overlap is exposed.
		TS_ASSERT( invdots[pro17_cb.rsd()]->atom_overlap_is_exposed( pro17_cb.atomno(), pro17_cd.atomno() ));
		//invdots[pro17_cb.rsd()]->write_circle_intersection_mask_to_kinemage( std::cout, pro17_cb.atomno(), *(invdots[ pro17_cd.rsd() ]), pro17_cd.atomno(), true );

		// the overlap between this atom pair is definitely buried. make sure buried overlap is called buried.
		TS_ASSERT( ! invdots[pro17_ca.rsd()]->atom_overlap_is_exposed( pro17_ca.atomno(), pro17_cd.atomno() ));
		//invdots[pro17_ca.rsd()]->write_circle_intersection_mask_to_kinemage( std::cout, pro17_ca.atomno(), *(invdots[ pro17_cd.rsd() ]), pro17_cd.atomno(), true );

		//TS_ASSERT( ! invdots[ile4_cg2.rsd()]->atom_overlap_is_exposed( ile4_cg2.atomno(), * invdots[8], lys8_cb.atomno() ));
	}


	void test_invrotamerdots_overlap_is_buried_1FKB() {
		TR << "Running test_invrotamerdots_overlap_is_buried_1FKB..." << std::endl;
		using namespace core::id;

		pose::Pose second_pose;
		core::import_pose::pose_from_file( second_pose, "core/pack/1FKB.pdb.gz" , core::import_pose::PDB_file);

		utility::vector1< RotamerDotsOP > rdots( 107 );
		utility::vector1< InvRotamerDotsOP > invdots( 107 );
		for ( Size ii = 1; ii <= 107; ++ii ) {
			rdots[ ii ] = RotamerDotsOP( new RotamerDots( second_pose.residue(ii).clone(), exclude_hydrogens, use_expanded_polar_atom_radii ) );
			rdots[ ii ]->increment_self_overlap();
		}
		for ( Size ii = 1; ii <= 107; ++ii ) {
			for ( Size jj = ii+1; jj <= 107; ++jj ) {
				rdots[ii]->increment_both( *rdots[jj] );
			}
			invdots[ii] = InvRotamerDotsOP( new InvRotamerDots() );
			invdots[ii]->setup_from_rotamer_dots( *rdots[ii] );
		}
		AtomID ala84_cb( second_pose.residue_type(84).atom_index( "CB" ), 84 );
		AtomID thr85_cg2( second_pose.residue_type(85).atom_index( "CG2" ), 85 );

		//invdots[84]->write_exposed_dots_to_kinemage( std::cout, true );
		//invdots[85]->write_exposed_dots_to_kinemage( std::cout, true );

		// the overlap between this atom pair is definitely buried. make sure buried overlap is called buried.
		TS_ASSERT( ! invdots[ ala84_cb.rsd() ]->atom_overlap_is_exposed( ala84_cb.atomno(), *( invdots[ thr85_cg2.rsd()] ), thr85_cg2.atomno() ));
		//invdots[ ala84_cb.rsd() ]->write_circle_intersection_mask_to_kinemage( std::cout, ala84_cb.atomno(), *(invdots[ thr85_cg2.rsd() ]), thr85_cg2.atomno(), true );

	}

};

