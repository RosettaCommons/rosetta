// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/denovo_design/BridgeChainsMover.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::connection::BridgeChainsMover
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/movers/BridgeChainsMover.hh>

// Protocol headers
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
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
#include <boost/assign.hpp>

// unit test utility functions
#include <protocols/denovo_design/test_utils.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.BridgeChainsMover.cxxtest" );

using namespace protocols::denovo_design::architects;
using namespace protocols::denovo_design::connection;
using namespace protocols::denovo_design::movers;

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
class BridgeChainsMoverTests : public CxxTest::TestSuite {
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
	void test_twohelix()
	{
		using namespace protocols::denovo_design;
		StructureDataFactory const & factory = *StructureDataFactory::get_instance();
		core::scoring::ScoreFunctionOP scorefxn = create_scorefxn();
		TS_ASSERT( scorefxn );

		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/denovo_design/connection/twohelix_structuredata.pdb" );
		core::pose::Pose const input_pose = pose;

		std::string const id_tag = "BridgeChainsMoverUnit";
		StructureData const orig = factory.get_from_pose( pose, id_tag );

		TS_ASSERT( factory.observer_attached( pose ) );
		TS_ASSERT_THROWS_NOTHING( orig.check_pose_consistency( pose ) );

		std::string const conn_id = "helixconn";
		std::string const s1_id = "BridgeChainsMoverUnit.H02";
		std::string const s2_id = "BridgeChainsMoverUnit.H01";
		BridgeChainsMover conn;
		conn.set_id( "helixconn" );
		conn.set_dry_run( true );
		conn.set_scorefxn( *scorefxn );
		conn.set_segment1_ids( s1_id );
		conn.set_segment2_ids( s2_id );
		conn.set_motifs( "2LX-1LA-1LG-1LB", "1:4" );
		conn.set_overlap( 0 );
		conn.apply( pose );

		StructureData perm = factory.get_from_const_pose( pose );
		TS_ASSERT_THROWS_NOTHING( perm.check_pose_consistency( pose ) );

		// segment1-->helixconn->segment2
		TS_ASSERT_EQUALS( input_pose.fold_tree().num_jump(), 1 );
		TS_ASSERT_EQUALS( input_pose.conformation().num_chains(), 2 );
		TS_ASSERT_EQUALS( pose.fold_tree().num_jump(), 0 );
		TS_ASSERT_EQUALS( pose.conformation().num_chains(), 1 );

		// check to make sure the right stuff is set
		TS_ASSERT_EQUALS( perm.segment( conn_id ).length(), 5 );
		TS_ASSERT_EQUALS( perm.segment( conn_id ).elem_length(), 5 );
		TS_ASSERT_EQUALS( perm.segment( conn_id ).cutpoint(), 0 );
		TS_ASSERT_EQUALS( perm.segment( conn_id ).ss(), "LLLLL" );
		TS_ASSERT_EQUALS( perm.segment( conn_id ).abego(), "XXAGB" );
		TS_ASSERT_EQUALS( perm.segment( s1_id ).upper_segment(), conn_id );
		TS_ASSERT_EQUALS( perm.segment( conn_id ).lower_segment(), s1_id );
		TS_ASSERT_EQUALS( perm.segment( conn_id ).upper_segment(), s2_id );
		TS_ASSERT_EQUALS( perm.segment( s2_id ).lower_segment(), conn_id );

		// check for unwanted movement
		check_unwanted_movement( orig, input_pose, perm, pose );
	}

	// test connectchains - interaction with components
	void test_conn_chains()
	{
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		core::pose::Pose input_pose;
		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/test_pdbcomp_BridgeChains.pdb" );
		core::pose::Pose const orig_pose = input_pose;
		StructureData const orig = StructureDataFactory::get_instance()->get_from_pose( input_pose );
		TS_ASSERT_THROWS_NOTHING( orig.check_pose_consistency( input_pose ) );
		TS_ASSERT_EQUALS( input_pose.conformation().num_chains(), 5 );

		// look for proper termini
		core::Size c_count = 0, n_count = 0;
		for ( core::Size i=1; i<=input_pose.total_residue(); ++i ) {
			if ( input_pose.residue(i).is_lower_terminus() ) {
				n_count += 1;
			}
			if ( input_pose.residue(i).is_upper_terminus() ) {
				c_count += 1;
			}
		}
		TS_ASSERT_EQUALS( n_count, 5 );
		TS_ASSERT_EQUALS( c_count, 5 );

		// create connectchains for testing purposes
		PoseArchitect pose_arch( "pose" );
		StructureDataOP perm = pose_arch.apply( input_pose );
		BridgeChainsMover conn;
		conn.set_id( "test" );
		conn.set_segment1_ids( "sheet1.s2" );
		conn.set_segment2_ids( "rot_comp.rot_comp2.catalytic.2" );
		conn.set_overlap( 1 );

		core::scoring::ScoreFunctionOP scorefxn = create_scorefxn();
		TS_ASSERT( scorefxn );
		conn.set_scorefxn( *scorefxn );
		conn.set_motifs( "1LE", "1" );
		core::Real random = 0.2501;
		TS_ASSERT_THROWS( conn.architect().apply( *perm->clone(), random ), EXCN_ConnectionSetupFailed );

		// now set it up with a longer loop that will build
		conn.set_motifs( "4LX", "3" );

		// do connecting
		conn.apply( input_pose );

		TS_ASSERT_EQUALS( conn.get_last_move_status(), protocols::moves::MS_SUCCESS );

		StructureData sd_conn = StructureDataFactory::get_instance()->get_from_pose( input_pose );
		sd_conn.check_pose_consistency( input_pose );

		for ( core::Size resid=1; resid<=input_pose.total_residue(); ++resid ) {
			TR.Debug << "After bridging: " << input_pose.residue(resid).name() << resid << std::endl;
		}
		TS_ASSERT_EQUALS( input_pose.conformation().num_chains(), 4 );
		TS_ASSERT_EQUALS( input_pose.fold_tree().num_jump(), 3 );

		// look for proper termini
		TR << sd_conn << std::endl;
		std::set< core::Size > const expected_lower_termini =
			boost::assign::list_of (1) (15) (35) (43);
		std::set< core::Size > const expected_upper_termini =
			boost::assign::list_of (14) (34) (42) (49);

		c_count = 0;
		n_count = 0;
		for ( core::Size resid=1; resid<=input_pose.total_residue(); ++resid ) {
			if ( input_pose.residue(resid).is_lower_terminus() ) {
				n_count += 1;
				TS_ASSERT_DIFFERS( expected_lower_termini.find( resid ), expected_lower_termini.end() );
			} else {
				TS_ASSERT_EQUALS( expected_lower_termini.find( resid ), expected_lower_termini.end() );
			}
			if ( input_pose.residue(resid).is_upper_terminus() ) {
				c_count += 1;
				TS_ASSERT_DIFFERS( expected_upper_termini.find( resid ), expected_upper_termini.end() );
			} else {
				TS_ASSERT_EQUALS( expected_upper_termini.find( resid ), expected_upper_termini.end() );
			}
		}
		TS_ASSERT_EQUALS( n_count, 4 );
		TS_ASSERT_EQUALS( c_count, 4 );

		// compare coordinates to make sure nothing changed that shouldn't have
		core::Size const ss = orig.lower_anchor("sheet1.s2");
		core::Size const se = orig.upper_anchor("sheet1.s2");
		core::Size const ms = orig.lower_anchor("rot_comp.rot_comp2.catalytic.2");
		core::Size const me = orig.upper_anchor("rot_comp.rot_comp2.catalytic.2");
		core::Vector xyz1 = orig_pose.residue(2).xyz("CA");
		core::Vector xyz2 = input_pose.residue( denovo_design::components::new_resnum( 2, ss, se, ms, me, 4 ) ).xyz("CA");
		for ( core::Size i=1; i<=orig_pose.total_residue(); ++i ) {
			core::Size const new_res( denovo_design::components::new_resnum( i, ss, se, ms, me, 4 ) );
			TR.Debug << "i=" << i << " new=" << new_res << std::endl;
			if ( ! new_res ) {
				continue;
			}

			// orig residues 22 and 40 are on the edges of remodeled area and can have small positional changes due to overlap == 1
			if ( i == 22 ) continue;
			if ( i == 40 ) continue;
			core::Vector const & orig_xyz = orig_pose.residue(i).xyz("CA");
			core::Vector const & newpt = input_pose.residue(new_res).xyz("CA");
			TS_ASSERT_DELTA( orig_xyz.distance(xyz1), newpt.distance(xyz2), 1e-6 );
			TS_ASSERT_DELTA( orig_pose.phi(i), input_pose.phi(new_res), 1e-6 );
			TS_ASSERT_DELTA( orig_pose.psi(i), input_pose.psi(new_res), 1e-6 );
		}
		TS_ASSERT( input_pose.fold_tree().check_fold_tree() );
	}


};
