// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/denovo_design/BridgeChains.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::connection::BridgeChains
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/connection/BridgeChains.hh>
#include <protocols/denovo_design/connection/Connection.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/util.hh>

// Protocol headers
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/ABEGOManager.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <numeric/random/random.hh>

// Boost headers
#include <boost/algorithm/string.hpp>

// unit test utility functions
#include <protocols/denovo_design/test_utils.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.connection.BridgeChains.cxxtest" );

namespace denovo_design {
namespace components {

/// @brief computes the new residue number after one component is joined to another using chain termini instead of anchors
core::Size new_resnum_termini(
	core::Size const orig_resnum,
	core::Size const css,
	core::Size const cse,
	core::Size const cms,
	core::Size const cme,
	core::Size const insert_length,
	core::Size const c1_cut_residues,
	core::Size const c2_cut_residues )
{
	if ( ( c2_cut_residues && ( orig_resnum == cms ) ) || ( c1_cut_residues && ( orig_resnum == cse ) ) ) {
		return 0;
	}
	// before movable segment
	if ( orig_resnum < cms ) {
		// before or inside fixed segment
		if ( orig_resnum <= cse ) {
			return orig_resnum;
			// after fixed segment
		} else {
			return orig_resnum + cme - cms - c2_cut_residues + insert_length;
		}
		// after movable segment
	} else if ( orig_resnum > cme ) {
		// before or inside fixed segment
		if ( orig_resnum <= cse ) {
			return orig_resnum - cme - 1 + cms;
			// after fixed segment
		} else {
			return orig_resnum - c1_cut_residues - c2_cut_residues + insert_length;
		}
		// inside movable segment
	} else {
		// before fixed segment
		if ( orig_resnum < css ) {
			return cse - cme - c1_cut_residues - c2_cut_residues + orig_resnum + insert_length;
			// after fixed segment
		} else if ( orig_resnum > cse ) {
			return cse - c1_cut_residues + orig_resnum - cms + insert_length;
			// also inside fixed segment -- should never happen
		} else {
			utility_exit_with_message( boost::lexical_cast< std::string >(orig_resnum) + " is inside both the fixed and movable segments. Something is wrong..." );
		}
		return cse - c1_cut_residues + orig_resnum - cms + insert_length;
	}
}

/// @brief computes the new residue number after one component is joined to another
core::Size new_resnum(
	core::Size const orig_resnum,
	core::Size const ss,
	core::Size const se,
	core::Size const ms,
	core::Size const me,
	core::Size const insert_length )
{
	return new_resnum_termini( orig_resnum, ss-1, se+1, ms-1, me+1, insert_length, 1, 1 );
}

}
}

// --------------- Test Class --------------- //
class BridgeChainsTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		protocols_init();

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	core::scoring::ScoreFunctionOP create_scorefxn() const
	{
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction() );
		scorefxn->set_weight( core::scoring::vdw, 1.0 );
		scorefxn->set_weight( core::scoring::rg, 1.0 );
		scorefxn->set_weight( core::scoring::rama, 0.1 );
		scorefxn->set_weight( core::scoring::hs_pair, 1.0 );
		scorefxn->set_weight( core::scoring::ss_pair, 1.0 );
		scorefxn->set_weight( core::scoring::rsigma, 1.0 );
		scorefxn->set_weight( core::scoring::coordinate_constraint, 1.0 );
		return scorefxn;
	}

	// test connectchains - interaction with components
	void test_conn_chains()
	{
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		using namespace protocols::denovo_design::connection;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/test_pdbcomp_BridgeChains.pdb" );
		StructureDataOP perm = StructureData::create_from_pose( input_pose, "BridgeChainsUnit" );
		TS_ASSERT( perm->pose() );
		TS_ASSERT_THROWS_NOTHING( perm->check_consistency() );
		TS_ASSERT_EQUALS( perm->pose()->conformation().num_chains(), 5 );

		// save original for testing
		StructureDataOP orig = perm->clone();

		// look for proper termini
		core::Size c_count = 0, n_count = 0;
		for ( core::Size i=1; i<=perm->pose()->total_residue(); ++i ) {
			if ( perm->pose()->residue(i).is_lower_terminus() ) {
				n_count += 1;
			}
			if ( perm->pose()->residue(i).is_upper_terminus() ) {
				c_count += 1;
			}
		}
		TS_ASSERT_EQUALS( n_count, 5 );
		TS_ASSERT_EQUALS( c_count, 5 );

		// create connectchains for testing purposes
		BridgeChains conn;
		conn.set_id( "test" );
		conn.set_comp1_ids( "sheet1.s1,sheet1.s2,sheet1.s3,sheet1.s4" );
		conn.set_comp2_ids( "rot_comp.rot_comp2.catalytic.2" );
		conn.set_overlap( 1 );
		conn.set_check_abego( false );

		core::scoring::ScoreFunctionOP scorefxn = create_scorefxn();
		TS_ASSERT( scorefxn );
		conn.set_scorefxn( scorefxn );
		conn.set_lengths( "1" );
		conn.set_cut_resis( "1" );
		try {
			// this should have failed since nothing is close enough
			conn.setup_from_random( *(perm->clone()), 0.2501 );
			TS_ASSERT( false );
		} catch ( EXCN_Setup const & e ) {
			TS_ASSERT( true );
		}

// now set it up with a longer loop that will build
		conn.set_lengths( "4" );
		conn.set_cut_resis( "3" );
		TR << "Setting up " << *perm << std::endl;
		conn.setup_from_random( *perm, 0.2501 );
		TS_ASSERT( conn.segments_fixed( *perm ) );
		TS_ASSERT( perm->pose() );

		// perform setup tasks (e.g. loop building)
		conn.setup( *perm );
		TS_ASSERT_EQUALS( conn.get_last_move_status(), protocols::moves::MS_SUCCESS );
		TS_ASSERT_EQUALS( perm->find_jump( conn.lower_segment_id(*perm) ), 0 );
		TS_ASSERT_EQUALS(
			perm->pose()->fold_tree().downstream_jump_residue( perm->find_jump( conn.upper_segment_id(*perm) ) ),
			perm->segment("rot_comp.rot_comp2.catalytic.2").safe() );

		// do connecting
		conn.apply_connection( *perm );
		TS_ASSERT_EQUALS( conn.get_last_move_status(), protocols::moves::MS_SUCCESS );

		// tell permutation we're connecting these termini
		perm->connect_segments( "test", "test_1" );
		perm->merge_segments( "test", "test_1", "test" );

		utility::vector1< std::string > roots;
		roots.push_back( conn.lower_segment_id(*perm) );
		roots.push_back( conn.upper_segment_id(*perm) );
		perm->consolidate_movable_groups( roots );
		TS_ASSERT_THROWS_NOTHING( perm->check_consistency() );
		for ( core::Size i=1; i<=perm->pose()->total_residue(); ++i ) {
			TR << perm->pose()->residue(i).name() << i << std::endl;
		}
		TS_ASSERT_EQUALS( perm->pose()->conformation().num_chains(), 4 );
		TS_ASSERT_EQUALS( perm->pose()->num_jump(), 3 );

		// look for proper termini
		c_count = 0;
		n_count = 0;
		for ( core::Size i=1; i<=perm->pose()->total_residue(); ++i ) {
			if ( perm->pose()->residue(i).is_lower_terminus() ) {
				n_count += 1;
			}
			if ( perm->pose()->residue(i).is_upper_terminus() ) {
				c_count += 1;
			}
		}
		TS_ASSERT_EQUALS( n_count, 4 );
		TS_ASSERT_EQUALS( c_count, 4 );

		core::pose::PoseOP new_pose = perm->pose()->clone();
		TS_ASSERT_EQUALS( perm->pose()->num_jump(), 3 );

		// look for proper termini
		c_count = 0;
		n_count = 0;
		for ( core::Size i=1; i<=perm->pose()->total_residue(); ++i ) {
			if ( perm->pose()->residue(i).is_lower_terminus() ) {
				n_count += 1;
			}
			if ( perm->pose()->residue(i).is_upper_terminus() ) {
				c_count += 1;
			}
		}
		TS_ASSERT_EQUALS( n_count, 4 );
		TS_ASSERT_EQUALS( c_count, 4 );

		// compare coordinates to make sure nothing changed that shouldn't have
		core::Size const ss = orig->lower_anchor("sheet1.s2");
		core::Size const se = orig->upper_anchor("sheet1.s2");
		core::Size const ms = orig->lower_anchor("rot_comp.rot_comp2.catalytic.2");
		core::Size const me = orig->upper_anchor("rot_comp.rot_comp2.catalytic.2");
		core::Vector xyz1 = orig->pose()->residue(2).xyz("CA");
		core::Vector xyz2 = new_pose->residue( denovo_design::components::new_resnum( 2, ss, se, ms, me, 4 ) ).xyz("CA");
		for ( core::Size i=1; i<=orig->pose()->total_residue(); ++i ) {
			core::Size const new_res( denovo_design::components::new_resnum( i, ss, se, ms, me, 4 ) );
			TR.Debug << "i=" << i << " new=" << new_res << std::endl;
			if ( ! new_res ) {
				continue;
			}
			// orig residues 22 and 40 are on the edges of remodeled area and can have small positional changes due to overlap == 1
			if ( ( i == 22 ) || ( i == 40 ) ) {
				continue;
			}
			core::Vector const & orig_xyz = orig->pose()->residue(i).xyz("CA");
			core::Vector const & newpt = perm->pose()->residue(new_res).xyz("CA");
			TS_ASSERT_DELTA( orig_xyz.distance(xyz1), newpt.distance(xyz2), 1e-6 );
			TS_ASSERT_DELTA( orig->pose()->phi(i), perm->pose()->phi(new_res), 1e-6 );
			TS_ASSERT_DELTA( orig->pose()->psi(i), perm->pose()->psi(new_res), 1e-6 );
		}
		TS_ASSERT( perm->pose()->fold_tree().check_fold_tree() );
	}

	// test connectchains - interaction with components
	void test_twohelix()
	{
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::connection;
		core::scoring::ScoreFunctionOP scorefxn = create_scorefxn();
		TS_ASSERT( scorefxn );

		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/twohelix_structuredata.pdb" );

		std::string const id_tag = "BridgeChainsUnit";
		components::StructureDataOP perm = components::StructureData::create_from_pose( input_pose, id_tag );
		TS_ASSERT( perm );
		TS_ASSERT( perm->pose() );
		TS_ASSERT_THROWS_NOTHING( perm->check_consistency() );

		// save original for comparison
		components::StructureDataOP orig = perm->clone();

		BridgeChains conn;
		conn.set_id( "helixconn" );
		conn.set_scorefxn( scorefxn );
		conn.set_comp1_ids( "2" );
		conn.set_comp2_ids( "1" );
		conn.set_lengths( "5" );
		conn.set_overlap( 2 );
		conn.setup_permutation( *perm );

		TS_ASSERT_THROWS_NOTHING( perm->check_consistency() );

		// segment1-->segment2, segment1-->helixconn, segment1-->helixconn_1
		TS_ASSERT_EQUALS( perm->pose()->fold_tree().num_jump(), 3 );
		conn.setup( *perm );

		// check to make sure the right stuff is set
		TS_ASSERT_EQUALS( conn.build_len(*perm), 5 );
		TS_ASSERT_EQUALS( perm->find_jump("2"), 0 );
		TS_ASSERT_EQUALS( perm->pose()->fold_tree().downstream_jump_residue( perm->find_jump("1") ),
			perm->segment( "1" ).safe() );
		TS_ASSERT_EQUALS( conn.build_abego(*perm), "XXXXX" );
		TS_ASSERT_EQUALS( conn.build_ss(*perm), "LLLLL" );
		TS_ASSERT_EQUALS( conn.build_len(*perm), 5 );
		TS_ASSERT_EQUALS( conn.user_chain1(), 0 );
		TS_ASSERT_EQUALS( conn.user_chain2(), 0 );
		TS_ASSERT_EQUALS( conn.upper_segment_id(*perm), "1" );
		TS_ASSERT_EQUALS( conn.lower_segment_id(*perm), "2" );

		// DSSP the shit out of this pose
		protocols::moves::DsspMover dssp;
		perm->apply_mover( dssp );

		// setup starting, ending residues
		core::Size start1, end1, start2, end2;
		int jump = 0;
		bool const_fold_tree( false );
		start1 = perm->segment(conn.lower_segment_id(*perm)).nterm_resi();
		end1 = perm->segment(conn.lower_segment_id(*perm)).cterm_resi();
		start2 = perm->segment(conn.upper_segment_id(*perm)).nterm_resi();
		end2 = perm->segment(conn.upper_segment_id(*perm)).cterm_resi();
		const_fold_tree = true;
		jump = perm->find_jump(conn.upper_segment_id(*perm));

		// delete unneeded terminal residues
		perm->delete_trailing_residues(conn.loop_lower(*perm));
		--start2;
		--end2;
		perm->delete_leading_residues(conn.loop_upper(*perm));
		--start2;
		--end2;

		TR << *perm << std::endl;
		TS_ASSERT_EQUALS( end1+1, perm->segment(conn.loop_lower(*perm)).nterm_resi() );
		TS_ASSERT_EQUALS( start2-1, perm->segment(conn.loop_upper(*perm)).cterm_resi() );
		TS_ASSERT_EQUALS( end1 + conn.build_len(*perm) + 1, start2 );

		TS_ASSERT_EQUALS( start1, perm->segment(conn.lower_segment_id(*perm)).nterm_resi() );
		TS_ASSERT_EQUALS( end1, perm->segment(conn.lower_segment_id(*perm)).cterm_resi() );
		TS_ASSERT_EQUALS( start2, perm->segment(conn.upper_segment_id(*perm)).nterm_resi() );
		TS_ASSERT_EQUALS( end2, perm->segment(conn.upper_segment_id(*perm)).cterm_resi() );
		TS_ASSERT( const_fold_tree );
		TS_ASSERT_EQUALS( jump, 1 );

		int const len( end2 - start2 + 1 );
		TS_ASSERT_EQUALS( len, 34 );

		// find left and right boundaries
		core::Size left( end1 - conn.lower_overlap() + 1 );
		if ( conn.lower_overlap() > end1 ) {
			left = 0;
		}
		core::Size right( start2 + conn.upper_overlap() - 1 );
		if ( left < 1 ) {
			left = 1;
		} else if ( left > perm->pose()->total_residue() ) {
			left = perm->pose()->total_residue();
		}
		if ( right < 1 ) {
			right = 1;
		} else if ( right > perm->pose()->total_residue() ) {
			right = perm->pose()->total_residue();
		}

		TS_ASSERT_EQUALS( perm->pose()->total_residue(), 73 );
		TS_ASSERT_EQUALS( left, 33 );
		TS_ASSERT_EQUALS( right, 41 );
		TS_ASSERT_EQUALS( perm->pose()->fold_tree().num_jump(), 1 );

		// test motif/parsing
		Connection::Motif motif = conn.parse_motif( "5:LX-3EB" );
		TS_ASSERT_EQUALS( motif.len, 8 );
		for ( core::Size i=1; i<=5; ++i ) {
			TS_ASSERT_EQUALS( motif.ss[i-1], 'L' );
			TS_ASSERT_EQUALS( motif.abego[i-1], 'X' );
		}
		for ( core::Size i=6; i<=8; ++i ) {
			TS_ASSERT_EQUALS( motif.ss[i-1], 'E' );
			TS_ASSERT_EQUALS( motif.abego[i-1], 'B' );
		}

		core::pose::Pose const & pose = *(perm->pose());
		// test SS
		std::string const ss = conn.ss_insert( pose, motif.ss, left, right, end1+1, start2-1 );
		TS_ASSERT_EQUALS( ss, "HHLLLLLEEEHH" );

		// test abegos
		utility::vector1< std::string > const complete_abego = core::sequence::get_abego( pose, 1 );
		utility::vector1< std::string > const abego = conn.abego_insert( complete_abego, motif.abego, left, right, end1+1, start2-1 );
		TS_ASSERT_EQUALS( abego.size(), 12 );
		TS_ASSERT_EQUALS( abego[1], "A" );
		TS_ASSERT_EQUALS( abego[2], "A" );
		TS_ASSERT_EQUALS( abego[3], "X" );
		TS_ASSERT_EQUALS( abego[4], "X" );
		TS_ASSERT_EQUALS( abego[5], "X" );
		TS_ASSERT_EQUALS( abego[6], "X" );
		TS_ASSERT_EQUALS( abego[7], "X" );
		TS_ASSERT_EQUALS( abego[8], "B" );
		TS_ASSERT_EQUALS( abego[9], "B" );
		TS_ASSERT_EQUALS( abego[10], "B" );
		TS_ASSERT_EQUALS( abego[11], "A" );
		TS_ASSERT_EQUALS( abego[12], "A" );
		check_unwanted_movement( *orig, *perm );
	}

};
