// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/pose/util.cxxtest.hh
/// @brief  unit tests for core::pose::util.hh/cc functions
/// @author James Thompson
/// @author Steven Lewis (compare coordinates)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>


#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class PoseUtilTests : public CxxTest::TestSuite {


public: // setup


	typedef core::Size Size;
	typedef core::pose::PDBInfo PDBInfo;
	typedef core::pose::PDBInfoOP PDBInfoOP;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;


public: //setup


	PoseUtilTests() {};


	// Shared initialization.
	void setUp() {
		core_init();
	}


	// Shared finalization.
	void tearDown() {}

public: // re-used methods

	/// @brief return a one-chain Pose ( 9 res )
	PoseOP one_chain_pose() {
		PoseOP pose( new Pose() );
		core::pose::make_pose_from_sequence(
			*pose,
			"A[ALA:NtermProteinFull]CDEFGHIK[LYS:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )
		);

		for ( core::Size i = 1; i <= pose->n_residue(); ++i ) {
			pose->set_secstruct( i, 'L' );
		}

		return pose;
	}


	/// @brief return a two-chain Pose ( 9 res + 11 res )
	PoseOP two_chain_pose() {
		PoseOP pose( new Pose() );
		core::pose::make_pose_from_sequence(
			*pose,
			"A[ALA:NtermProteinFull]CDEFGHIK[LYS:CtermProteinFull]L[LEU:NtermProteinFull]MNPQRSTVWY[TYR:CtermProteinFull]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )
		);

		for ( core::Size i = 1; i <= 20; ++i ) {
			pose->set_secstruct( i, 'L' );
		}

		return pose;
	}


	PDBInfoOP add_empty_pdb_info( Pose & pose ) {
		PDBInfoOP pdbinfo( new PDBInfo( pose.n_residue() ) );
		pose.pdb_info( pdbinfo );

		return pdbinfo;
	}


public: // tests


	/// @brief test Pose DataCache manipulation methods
	void test_pose_string_map() {
		using namespace core::scoring;

		Pose pose;

		// test string-based comments
		std::string const key( "key" );
		std::string const desired_val( "not_empty" );
		std::string val( "empty" );
		bool key_exists;

		key_exists = get_comment( pose, key, val );
		TS_ASSERT( key_exists == false );
		TS_ASSERT( val == "empty" );

		add_comment( pose, key, desired_val );

		key_exists = get_comment( pose, key, val );
		TS_ASSERT( key_exists == true );
		TS_ASSERT( val == desired_val );

		std::map< std::string, std::string > comments = get_all_comments( pose );
		std::map< std::string, std::string >::const_iterator it
			= comments.find( key );

		TS_ASSERT( it != comments.end() );
		TS_ASSERT( it->first  == key );
		TS_ASSERT( it->second == val );
	} // test_pose_string_map


	void test_pose_float_map() {
		using namespace core::scoring;
		using namespace core::pose;

		Pose pose;

		// test string-based comments
		std::string const key( "key" );
		float const desired_val( 31882.0 );
		core::Real val( 0.0 );
		bool key_exists;

		key_exists = getPoseExtraScore( pose, key, val );
		TS_ASSERT( key_exists == false );
		TS_ASSERT( val == 0.0 );

		setPoseExtraScore( pose, key, desired_val );

		key_exists = getPoseExtraScore( pose, key, val );
		TS_ASSERT( key_exists == true );
		TS_ASSERT( val == desired_val );

		val = 0.0;
		clearPoseExtraScore( pose, key );
		key_exists = getPoseExtraScore( pose, key, val );
		TS_ASSERT( key_exists == false );
		TS_ASSERT( val == 0.0 )
	} // test_pose_float_map


	/// @brief test renumber_pdbinfo_based_on_conf_chains()
	void test_renumber_pdbinfo_based_on_conf_chains() {
		using core::pose::renumber_pdbinfo_based_on_conf_chains;

		PoseOP pose_one = one_chain_pose();
		PoseOP pose_two = two_chain_pose();

		// TEST: single pose case with fully empty records
		add_empty_pdb_info( *pose_one );
		renumber_pdbinfo_based_on_conf_chains( *pose_one, true, true, false, false );

		for ( Size i = 1, ie = pose_one->n_residue(); i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_one->pdb_info()->number( i ), static_cast< int >( i ) );
			TS_ASSERT_EQUALS( pose_one->pdb_info()->chain( i ), 'A' );
		}

		// TEST: keeping existing starting number
		add_empty_pdb_info( *pose_two );
		pose_two->pdb_info()->chain( 1 , 'A' );
		pose_two->pdb_info()->number( 1 , 27 );
		pose_two->pdb_info()->chain( 10, 'B' );
		pose_two->pdb_info()->number( 10, -4 );
		renumber_pdbinfo_based_on_conf_chains( *pose_two, true, true, false, false );

		for ( Size i = 1, ie = 9; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i + 26 ) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), 'A' );
		}

		for ( Size i = 10, ie = 20; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i - 9 - 5) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), 'B' );
		}

		// TEST: keeping icodes
		add_empty_pdb_info( *pose_two );
		pose_two->pdb_info()->icode( 4, 'A' );
		pose_two->pdb_info()->icode( 5, 'B' );
		renumber_pdbinfo_based_on_conf_chains( *pose_two, true, true, true, false );

		for ( Size i = 1, ie = 3; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i ) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), 'A' );
		}

		// next two should have pdb_seq = 3 due to insertion codes
		for ( Size i = 4, ie = 5; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), 3 );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), 'A' );
		}

		for ( Size i = 6, ie = 9; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i - 2 ) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), 'A' );
		}

		for ( Size i = 10, ie = 20; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i - 9 ) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), 'B' );
		}

		// TEST: resetting icodes
		add_empty_pdb_info( *pose_two );
		pose_two->pdb_info()->icode( 4, 'A' );
		pose_two->pdb_info()->icode( 5, 'B' );
		renumber_pdbinfo_based_on_conf_chains( *pose_two, true, true, false, false );

		TS_ASSERT_EQUALS( pose_two->pdb_info()->icode( 4 ), ' ' );
		TS_ASSERT_EQUALS( pose_two->pdb_info()->icode( 5 ), ' ' );
		for ( Size i = 1, ie = 9; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i ) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), 'A' );
		}

		for ( Size i = 10, ie = 20; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i - 9 ) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), 'B' );
		}

		// TEST: no fix chains + icode resetting
		add_empty_pdb_info( *pose_two );
		pose_two->pdb_info()->icode( 4, 'A' );
		pose_two->pdb_info()->icode( 5, 'B' );
		renumber_pdbinfo_based_on_conf_chains( *pose_two, false, true, false, false );

		for ( Size i = 1, ie = 9; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i ) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), PDBInfo::empty_record() );
		}

		for ( Size i = 10, ie = 20; i <= ie; ++i ) {
			TS_ASSERT_EQUALS( pose_two->pdb_info()->number( i ), static_cast< int >( i - 9 ) );
			TS_ASSERT_EQUALS( pose_two->pdb_info()->chain( i ), PDBInfo::empty_record() );
		}

	} // test_renumber_pdbinfo_based_on_conf_chains

	/// @brief test compare_atom_coordinates
	void test_compare_atom_coordinates() {

		//6 digits of comparison precision
		core::Size const precision(6);

		core::pose::Pose const p1 = create_twores_1ubq_pose();
		core::pose::Pose p2(p1);

		//is a pose equal to itself?
		TS_ASSERT(core::pose::compare_atom_coordinates(p1, p1, precision));
		TS_ASSERT(core::pose::compare_atom_coordinates(p2, p2, precision));

		//is a pose equal to its copy?
		TS_ASSERT(core::pose::compare_atom_coordinates(p2, p1, precision));
		TS_ASSERT(core::pose::compare_atom_coordinates(p1, p2, precision));

		//test failure conditions
		//fail because lengths are different (special case) - this is expected to produce a warning
		p2.delete_polymer_residue(2);
		TS_ASSERT(!core::pose::compare_atom_coordinates(p2, p1, precision));
		TS_ASSERT(!core::pose::compare_atom_coordinates(p1, p2, precision));
		p2 = p1;

		//test nonmatching ResidueTypes (special case) - this is expected to produce a warning
		p2.replace_residue(2, *(core::conformation::ResidueFactory::create_residue(core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("SER"))), true);
		TS_ASSERT(!core::pose::compare_atom_coordinates(p2, p1, precision));
		TS_ASSERT(!core::pose::compare_atom_coordinates(p1, p2, precision));
		p2 = p1;

		//test nonmatching numbers of atoms (special case)
		//I have no idea how to do this - it shouldn't be possible given the ResidueType check

		//test nonmatching atomic coordinates (normal case)
		p2.set_psi(1, p1.psi(1)+15);
		TS_ASSERT(!core::pose::compare_atom_coordinates(p2, p1, precision));
		TS_ASSERT(!core::pose::compare_atom_coordinates(p1, p2, precision));

	}

	void test_conformation_sha1() {
		core::pose::Pose const p1 = create_pdb_string_2res_1ten_2res_trp_cage_pose();
		core::pose::Pose p2(p1);

		std::string p1_hash(core::pose::get_sha1_hash_excluding_chain('A',p1));
		std::string p2_hash(core::pose::get_sha1_hash_excluding_chain('A',p2));
		TS_ASSERT(p1_hash == p2_hash);
	}

}; // class PoseUtilTests
