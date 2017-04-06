// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

#include <core/kinematics/FoldTree.hh>

#include <basic/Tracer.hh>

#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR("core.pose.util.cxxtest");


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

		for ( core::Size i = 1; i <= pose->size(); ++i ) {
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
		PDBInfoOP pdbinfo( new PDBInfo( pose.size() ) );
		pose.pdb_info( pdbinfo );

		return pdbinfo;
	}


public: // tests

	void test_fix_pdbinfo_damaged_by_insertion() {
		// True PDBInfo based test
		Pose pose = *one_chain_pose();
		add_empty_pdb_info( pose );
		// Set up PDBInfo as though residue 5 is an insertion
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( ii == 5 ) {
				pose.pdb_info()->number( ii, 0 );
				pose.pdb_info()->chain( ii, ' ' );
				pose.pdb_info()->icode( ii, ' ' );
				continue;
			}
			pose.pdb_info()->number( ii, ii );
			pose.pdb_info()->chain( ii, 'A' );
			pose.pdb_info()->icode( ii, ' ' );
		}
		fix_pdbinfo_damaged_by_insertion( pose );
		TS_ASSERT( pose.pdb_info()->number( 5 ) == 5  ) ;
		TS_ASSERT( pose.pdb_info()->chain( 5 )  == 'A' );
		TS_ASSERT( pose.pdb_info()->icode( 5 )  == ' ' );
		TS_ASSERT( pose.pdb_info()->number( 6 ) == 6  ) ;
		TS_ASSERT( pose.pdb_info()->chain( 6 )  == 'A' );
		TS_ASSERT( pose.pdb_info()->icode( 6 )  == ' ' );

		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( ii == 5 ) {
				pose.pdb_info()->number( ii, 0 );
				pose.pdb_info()->chain( ii, ' ' );
				pose.pdb_info()->icode( ii, ' ' );
			} else if ( ii > 5 ) {
				pose.pdb_info()->number( ii, ii - 1 );
				pose.pdb_info()->chain( ii, 'A' );
				pose.pdb_info()->icode( ii, ' ' );
			} else {
				pose.pdb_info()->number( ii, ii );
				pose.pdb_info()->chain( ii, 'A' );
				pose.pdb_info()->icode( ii, ' ' );
			}
		}
		fix_pdbinfo_damaged_by_insertion( pose );
		TS_ASSERT( pose.pdb_info()->number( 5 ) == 4  ) ;
		TS_ASSERT( pose.pdb_info()->chain( 5 )  == 'A' );
		TS_ASSERT( pose.pdb_info()->icode( 5 )  == 'A' );
		TS_ASSERT( pose.pdb_info()->number( 6 ) == 5  ) ;
		TS_ASSERT( pose.pdb_info()->chain( 6 )  == 'A' );
		TS_ASSERT( pose.pdb_info()->icode( 6 )  == ' ' );

		/*Pose pose = *two_chain_pose();
		add_empty_pdb_info( pose );
		renumber_pdbinfo_based_on_conf_chains( pose, true, true, false, false );
		ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
		ResidueCOP res = rts->name_map( "ALA" );
		Pose test_one_pose = pose;
		pose.conformation().append_polymer_residue_after_seqpos( *rsd, 5, true );

		TS_ASSERT( test_one_pose.pdb_info()->number( 6 ) == 0  ) ;
		TS_ASSERT( test_one_pose.pdb_info()->chain( 6 )  == ' ' );
		TS_ASSERT( test_one_pose.pdb_info()->icode( 6 )  == ' ' );
		fix_pdbinfo_damaged_by_insertion( pose );
		TS_ASSERT( test_one_pose.pdb_info()->number( 6 ) == 5  ) ;
		TS_ASSERT( test_one_pose.pdb_info()->chain( 6 )  == 'A' );
		TS_ASSERT( test_one_pose.pdb_info()->icode( 6 )  == 'A' );
		*/
	}


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
		TS_ASSERT( val == 0.0 );
	} // test_pose_float_map


	/// @brief test renumber_pdbinfo_based_on_conf_chains()
	void test_renumber_pdbinfo_based_on_conf_chains() {
		using core::pose::renumber_pdbinfo_based_on_conf_chains;

		PoseOP pose_one = one_chain_pose();
		PoseOP pose_two = two_chain_pose();

		// TEST: single pose case with fully empty records
		add_empty_pdb_info( *pose_one );
		renumber_pdbinfo_based_on_conf_chains( *pose_one, true, true, false, false );

		for ( Size i = 1, ie = pose_one->size(); i <= ie; ++i ) {
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

	void test_unique_chains() {

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pose/2WDQ__tr.pdb" , core::import_pose::PDB_file);

		utility::vector1< bool > uniq_chains( compute_unique_chains( pose ) );

		TS_ASSERT( uniq_chains[1] == false );
		TS_ASSERT( uniq_chains[588] == false );
		TS_ASSERT( uniq_chains[589] == false );
		TS_ASSERT( uniq_chains[826] == false );

		TS_ASSERT( uniq_chains[827] == false );
		TS_ASSERT( uniq_chains[947] == false );
		TS_ASSERT( uniq_chains[948] == false );
		TS_ASSERT( uniq_chains[1052] == false );

		TS_ASSERT( uniq_chains[1053] == false );
		TS_ASSERT( uniq_chains[1640] == false );
		TS_ASSERT( uniq_chains[1641] == false );
		TS_ASSERT( uniq_chains[1878] == false );

		TS_ASSERT( uniq_chains[1879] == false );
		TS_ASSERT( uniq_chains[1999] == false );
		TS_ASSERT( uniq_chains[2000] == false );
		TS_ASSERT( uniq_chains[2104] == false );

		TS_ASSERT( uniq_chains[2105] == true );
		TS_ASSERT( uniq_chains[2692] == true );
		TS_ASSERT( uniq_chains[2693] == true );
		TS_ASSERT( uniq_chains[2930] == true );

		TS_ASSERT( uniq_chains[2931] == true );
		TS_ASSERT( uniq_chains[3051] == true );
		TS_ASSERT( uniq_chains[3052] == true );
		TS_ASSERT( uniq_chains[3156] == true );
	}

	//////////////////////////////////////////
	// Chain/jump methods
	//
public: // Utility methods
	int get_parent_jump( core::pose::Pose const & pose, core::Size resi ) const {
		core::kinematics::FoldTree const & fold_tree( pose.fold_tree() );
		int curr_pos( resi ), parent( resi );
		bool connected_by_jump( false );
		while ( parent != 0 && ! connected_by_jump ) {
			curr_pos = parent;
			parent = fold_tree.get_parent_residue( curr_pos, connected_by_jump );
		}
		if ( parent == 0 ) {
			// get_parent_residue() return 0 for the root atom
			return 0; // jump 0 is no jump.
		}
		// parent is the parent of curr_pos, connected via a jump.
		core::Size jumpno( fold_tree.jump_nr(parent, curr_pos ) );
		debug_assert( jumpno != 0 );
		return jumpno;
	}

	core::pose::PoseOP get_mixedup_pose() const {
		PoseOP pose( new Pose() );
		core::pose::make_pose_from_sequence(
			*pose,
			"X[VRT]/A[ALA:NtermProteinFull]CD/Z[MG]/F[PHE:NtermProteinFull]GHIKL[LEU:CtermProteinFull]/Z[MG]/Z[MG]/NPR[ARG:CtermProteinFull]/X[VRT]/Z[MG]",
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);
		if ( ! pose->pdb_info() ) { pose->pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( *pose ) ) ); }
		pose->pdb_info()->set_chains( std::string( "XAAAMBBBBCBMMAAAXM" ) );
		// Strap in: pos 1 is a virtual root, which has jumps to chain A, B, and the virtual root for the metals
		// Chain A is discontinous, and is connected by a jump in the middle
		// Chain B is contingous, but has a jump in the middle (note that this doesn't reset the chain id)
		//    It also has a residue labeled as "chain C", but actually part of the same polymer (jump and chain ID) as the rest of chain B
		// The Mg are all chain M, but most are connected to a virtual as jumps
		// The execeptions is the two Mg after chain B, which are connected to each other in backwards order,
		// and the very last Mg, which is connected to the end of chain B, and on a completely separate branch of the fold tree
		core::kinematics::FoldTree ft(pose->size());
		ft.clear();
		ft.add_edge( 1,  2,  1);
		ft.add_edge( 2,  4,  core::kinematics::Edge::PEPTIDE);
		ft.add_edge( 4, 14,  9); // chain A internal
		ft.add_edge(14, 16,  core::kinematics::Edge::PEPTIDE);
		ft.add_edge( 1,  6,  2);
		ft.add_edge( 6,  8,  core::kinematics::Edge::PEPTIDE);
		ft.add_edge( 8,  9,  8); // chain B internal
		ft.add_edge( 9, 11,  core::kinematics::Edge::PEPTIDE);
		ft.add_edge( 1, 17,  3);
		ft.add_edge(11, 18,  4); // End of chain B to last MG
		ft.add_edge(17, 13,  5);
		ft.add_edge(17,  5,  6);
		ft.add_edge(13, 12,  7);
		debug_assert( ft.connected() );
		debug_assert( ft.check_fold_tree() );
		pose->fold_tree( ft );
		return pose;
	}

public: // tests

	void test_chain_designations() {
		PoseOP pose( get_mixedup_pose() );
		core::pose::PDBInfoOP pdb_info( pose->pdb_info() );
		TS_ASSERT( pdb_info );

		utility::vector1< char > chain_letters{ 'X', 'A', 'A', 'A', 'M', 'B', 'B', 'B', 'B', 'C', 'B', 'M', 'M', 'A', 'A', 'A', 'X', 'M' };
		debug_assert( pose->size() == chain_letters.size() );
		// Chain IDs are always sequential.
		utility::vector1< core::Size > chain_ids{ 1, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 5, 6, 7, 7, 7, 8, 9 };
		debug_assert( pose->size() == chain_ids.size() );
		utility::vector1< core::Size > jumps{ 0, 1, 1, 1, 6, 2, 2, 2, 8, 8, 8, 7, 5, 9, 9, 9, 3, 4 };
		debug_assert( pose->size() == chain_letters.size() );
		for ( core::Size ii(1); ii <= pose->size(); ++ii ) {
			TR << "Residue " << ii
				<< " letter: " << pdb_info->chain( ii )
				<< " chain id: " << pose->chain( ii )
				<< " jump: " << get_parent_jump( *pose, ii ) << std::endl;
			TS_ASSERT_EQUALS( pdb_info->chain( ii ), chain_letters[ii] );
			TS_ASSERT_EQUALS( pose->chain( ii ), chain_ids[ii] );
			TS_ASSERT_EQUALS( get_parent_jump( *pose, ii ), jumps[ii] );
		}
	}

	void test_jump_functions() {
		PoseOP pose( get_mixedup_pose() );

		// jumps_from_pose  -- Rather silly function - just get a set with numbers 1 to numjumps
		core::pose::Jumps jumps;
		core::pose::jumps_from_pose( *pose, jumps );
		TS_ASSERT_EQUALS( jumps.size(), 9 );
		TS_ASSERT_EQUALS( jumps, core::pose::Jumps({ 1, 2, 3, 4, 5, 6, 7, 8, 9 }) );

		// get_jump_id_from_chain_id -- the jump(s) *directly* upstream of the chain_id
		TS_ASSERT_EQUALS( get_jump_ids_from_chain_ids( std::set<core::Size>{2,5,7,4}, *pose ), (std::set<core::Size>{1,7,9,2,8}) );
		TS_ASSERT_EQUALS( get_jump_id_from_chain_id( 9, *pose ), 4 );
		TS_ASSERT_EQUALS( get_jump_id_from_chain_id( 7, *pose ), 9 );
		TS_ASSERT_EQUALS( get_jump_id_from_chain_id( 2, *pose ), 1 );

		// get_jump_ids_from_chain
		TS_ASSERT_EQUALS( get_jump_ids_from_chain( "A", *pose ), (utility::vector1<core::Size>{1,9}) );
		TS_ASSERT_EQUALS( get_jump_ids_from_chain( 'B', *pose ), (utility::vector1<core::Size>{2,8}) );
		// C is a hard exit
		TS_ASSERT_EQUALS( get_jump_ids_from_chain( "M", *pose ), (utility::vector1<core::Size>{4,5,6,7}) );
		TS_ASSERT_EQUALS( get_jump_ids_from_chain( 'X', *pose ).size(), 1 ); // Root doesn't get a jump
		TS_ASSERT_EQUALS( get_jump_ids_from_chain( 'X', *pose ), (utility::vector1<core::Size>{3}) ); // Root doesn't get a jump
		TS_ASSERT_EQUALS( get_jump_ids_from_chain( 'T', *pose ), (utility::vector1<core::Size>{}) );

		TS_ASSERT_EQUALS( get_jump_id_from_chain( 'A', *pose ), 1 );
		TS_ASSERT_EQUALS( get_jump_id_from_chain( "B", *pose ), 2 );
		// C is a hard exit
		TS_ASSERT_EQUALS( get_jump_id_from_chain( "X", *pose ), 3 );
		TS_ASSERT_EQUALS( get_jump_id_from_chain( 'M', *pose ), 4 );

	}

	void test_chainid_functions() {
		PoseOP pose( get_mixedup_pose() );

		// get_chains -- Should just be a vector with serial numbers of chains, more or less
		utility::vector1< int > chains( core::pose::get_chains(*pose) );
		std::sort( chains.begin(), chains.end() );
		TS_ASSERT_EQUALS( chains.size(), 9 );
		TS_ASSERT_EQUALS( chains, utility::vector1< int >({ 1, 2, 3, 4, 5, 6, 7, 8, 9 }) );

		// conf2pdb_chain -- Not strictly speaking good on the mixed-up pose, as chain 4 is not consistent (it has both B and C chains
		// This will flag a warning on the tracer and return a naive assignment.
		// To test properly, we cheat by making the pose consistent and then putting it back
		TS_ASSERT_EQUALS( pose->pdb_info()->chain( 10 ), 'C' );
		pose->pdb_info()->chain( 10, 'B' );
		std::map< core::Size, char > chainid_to_chanlett( core::pose::conf2pdb_chain(*pose) );
		TS_ASSERT_EQUALS( chainid_to_chanlett.size(), 9 );
		TS_ASSERT_EQUALS( chainid_to_chanlett, (std::map< core::Size, char >( { {1,'X'}, {2,'A'}, {3,'M'}, {4,'B'}, {5,'M'}, {6,'M'}, {7,'A'}, {8,'X'}, {9,'M'} } )) );
		pose->pdb_info()->chain( 10, 'C' ); // Put it back, properly mixed up.

		// chain_end_res
		utility::vector1< core::Size > expected_chain_ends{ 1, 4, 5, 11, 12, 13, 16, 17, 18 };
		utility::vector1< core::Size > chain_ends( core::pose::chain_end_res( *pose ) );
		TS_ASSERT_EQUALS( chain_ends.size(), expected_chain_ends.size() );
		TS_ASSERT_EQUALS( chain_ends, expected_chain_ends );
		for ( core::Size ii(1); ii <= expected_chain_ends.size(); ++ii ) {
			TS_ASSERT_EQUALS( core::pose::chain_end_res( *pose, ii ), expected_chain_ends[ ii ] );
		}

		// core::pose::compute_unique_chains  -- Tested above.

		// has_chain
		for ( core::Size ii(1); ii <= pose->conformation().num_chains(); ++ii ) {
			TS_ASSERT( core::pose::has_chain( ii, *pose ) );
		}
		TS_ASSERT( ! core::pose::has_chain( core::Size(12), *pose ) );
		TS_ASSERT( ! core::pose::has_chain( core::Size(100), *pose ) );

		// get_resnums_for_chain_id
		TS_ASSERT_EQUALS( get_resnums_for_chain_id( *pose,  2 ), (utility::vector1< core::Size >{2,3,4}) );
		TS_ASSERT_EQUALS( get_resnums_for_chain_id( *pose,  4 ), (utility::vector1< core::Size >{6,7,8,9,10,11}) );
		TS_ASSERT_EQUALS( get_resnums_for_chain_id( *pose,  7 ), (utility::vector1< core::Size >{14,15,16}) );
		// TS_ASSERT_EQUALS( get_resnums_for_chain_id( *pose, 12 ), (utility::vector1< core::Size >{}) ); // debug_assert error

		// get_chain_residues
		TS_ASSERT_EQUALS( get_chain_residues( *pose,  2 ).size(), 3 );
		TS_ASSERT_EQUALS( get_chain_residues( *pose,  4 ).size(), 6 );
		TS_ASSERT_EQUALS( get_chain_residues( *pose,  8 ).size(), 1 );
		// TS_ASSERT_EQUALS( get_chain_residues( *pose, 12 ).size(), 0 ); // debug_assert error

		// get_residues_from_chains
		TS_ASSERT_EQUALS( get_residues_from_chains( *pose, utility::vector1<core::Size>{2} ).size(), 3 );
		TS_ASSERT_EQUALS( get_residues_from_chains( *pose, utility::vector1<core::Size>{3,9,6,5} ).size(), 4 );
		TS_ASSERT_EQUALS( get_residues_from_chains( *pose, utility::vector1<core::Size>{4,2,7} ).size(), 12 );
		//TS_ASSERT_EQUALS( get_residues_from_chains( *pose, utility::vector1<core::Size>{12,100} ).size(), 0 ); //debug_assert

		// get_chain_ids_from_chains
		TS_ASSERT_EQUALS( get_chain_ids_from_chains( utility::vector1<std::string>{"A","M"}, *pose ), (utility::vector1<core::Size>{2,3,5,6,7,9}) );
		TS_ASSERT_EQUALS( get_chain_ids_from_chains( utility::vector1<char>{'B','X'}, *pose ), (utility::vector1<core::Size>{1,4,8}) );
		TS_ASSERT_EQUALS( get_chain_ids_from_chains( utility::vector1<char>{'C'}, *pose ), (utility::vector1<core::Size>{4}) );
		TS_ASSERT_EQUALS( get_chain_ids_from_chains( utility::vector1<std::string>{"T"}, *pose ), (utility::vector1<core::Size>{}) );

		TS_ASSERT_EQUALS( get_chain_ids_from_chain( 'A', *pose), (utility::vector1<core::Size>{2,7}) );
		TS_ASSERT_EQUALS( get_chain_ids_from_chain( 'B', *pose), (utility::vector1<core::Size>{4}) );
		TS_ASSERT_EQUALS( get_chain_ids_from_chain( "C", *pose), (utility::vector1<core::Size>{4}) );
		TS_ASSERT_EQUALS( get_chain_ids_from_chain( "M", *pose), (utility::vector1<core::Size>{3,5,6,9}) );
		TS_ASSERT_EQUALS( get_chain_ids_from_chain( 'T', *pose), (utility::vector1<core::Size>{}) );

		TS_ASSERT_EQUALS( get_chain_id_from_chain( 'B', *pose), 4);
		TS_ASSERT_EQUALS( get_chain_id_from_chain( "C", *pose), 4);

		// get_chain_id_from_jump_id
		utility::vector1< core::Size > jump_id_to_chain_id{ 2, 4, 8, 9, 6, 3, 5, 4, 7 };
		debug_assert( jump_id_to_chain_id.size() == pose->num_jump() );
		for ( core::Size ii(1); ii <= pose->num_jump(); ++ii ) {
			TS_ASSERT_EQUALS( get_chain_id_from_jump_id( ii, *pose ), jump_id_to_chain_id[ii] );
		}
	}

	void test_chainletter_functions() {
		PoseOP pose( get_mixedup_pose() );

		// has_chain
		TS_ASSERT( core::pose::has_chain( 'A', *pose ) );
		TS_ASSERT( core::pose::has_chain( 'B', *pose ) );
		TS_ASSERT( core::pose::has_chain( 'C', *pose ) );
		TS_ASSERT( ! core::pose::has_chain( 'D', *pose ) );
		TS_ASSERT( core::pose::has_chain( 'M', *pose ) );
		TS_ASSERT( core::pose::has_chain( 'X', *pose ) );
		TS_ASSERT( ! core::pose::has_chain( 'Z', *pose ) );

		// res_in_chain
		utility::vector1< std::string > chain_letters{ "X", "A", "A", "A", "M", "B", "B", "B", "B", "C", "B", "M", "M", "A", "A", "A", "X", "M" };
		utility::vector1< std::string >    not_chains{ "A", "B", "M", "C", "A", "A", "C", "X", "M", "B", "C", "B", "C", "C", "M", "H", "H", "H" };
		TS_ASSERT_EQUALS( pose->size(), chain_letters.size() );
		TS_ASSERT_EQUALS( pose->size(), not_chains.size() );
		for ( core::Size ii(1); ii <= pose->size(); ++ii ) {
			TS_ASSERT( core::pose::res_in_chain( *pose, ii, chain_letters[ii] ) );
			TS_ASSERT( ! core::pose::res_in_chain( *pose, ii, not_chains[ii] ) );
		}

		// get_resnums_for_chain
		TS_ASSERT_EQUALS( get_resnums_for_chain( *pose, 'A' ), (utility::vector1< core::Size >{2,3,4,14,15,16}) );
		TS_ASSERT_EQUALS( get_resnums_for_chain( *pose, 'M' ), (utility::vector1< core::Size >{5,12,13,18}) );
		TS_ASSERT_EQUALS( get_resnums_for_chain( *pose, 'C' ), (utility::vector1< core::Size >{10}) );
		TS_ASSERT_EQUALS( get_resnums_for_chain( *pose, 'D' ), (utility::vector1< core::Size >{}) );

		// get_chain_from_chain_id -- this is the chain letter of the first residue in the chain number
		utility::vector1<char> ids_to_letters{ 'X', 'A', 'M', 'B', 'M', 'M', 'A', 'X', 'M' };
		debug_assert( ids_to_letters.size() == pose->num_chains() );
		for ( core::Size ii(1); ii <= pose->num_chains(); ++ii ) {
			TS_ASSERT_EQUALS( get_chain_from_chain_id( ii, *pose ), ids_to_letters[ii] );
		}

		// get_chain_from_jump_id
		utility::vector1< char > jump_id_to_chain{ 'A', 'B', 'X', 'M', 'M', 'M', 'M', 'B', 'A' };
		debug_assert( jump_id_to_chain.size() == pose->num_jump() );
		for ( core::Size ii(1); ii <= pose->num_jump(); ++ii ) {
			TS_ASSERT_EQUALS( get_chain_from_jump_id( ii, *pose ), jump_id_to_chain[ii] );
		}
	}

}; // class PoseUtilTests

