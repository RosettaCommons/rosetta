// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/CCDLoopClosureMover.cxxtest.hh
/// @brief test suite for protocols/loops/loop_closure/ccd/CCDLoopClosureMover
/// @author Brian D. Weitzner
/// @author Jason W. Labonte

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/RamaCheck.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/moves/MoverFactory.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// For RosettaScripts testing
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>

// Utility header
#include <utility/excn/EXCN_Base.hh>
#include <utility/string_util.hh>

// C++ header
#include <iostream>

namespace {

using core::id::BB;
using core::id::TorsionID;
using core::pose::Pose;
using core::pose::PoseOP;
using protocols::loops::Loop;
using core::kinematics::MoveMap;
using core::kinematics::MoveMapOP;
using protocols::loops::loop_closure::ccd::CCDLoopClosureMover;

static THREAD_LOCAL basic::Tracer TR( "test.protocols.loops.closers.CCDLoopClosureMover.cxxtest" );

class TestCCDLoopClosureMover : public CxxTest::TestSuite {

private:
	PoseOP pose_;

public:
	void setUp() {
		protocols_init();
		pose_ = core::import_pose::pose_from_file( "protocols/loops/2GB3.pdb" , core::import_pose::PDB_file);
	}

	// TODO: test both flavors of rama check -- this should probably be a separate unit test suite.

	// give an example loop and close it -- make sure the max cycles aren't hit, make sure the deviation is acceptable
	void test_CCDLoopClosureMover_closes() {
		Pose pose( *pose_ );

		// make sure we don't move omega
		MoveMapOP mm( new MoveMap );
		mm->set_bb( true );
		for ( core::uint i = 1; i < pose.total_residue(); ++i ) {
			mm->set( TorsionID( i, BB, 3 ), false );
		}

		Loop loop( 7, 14, 10 );

		set_single_loop_fold_tree( pose, loop );
		add_single_cutpoint_variant( pose, loop );

		for ( core::uint i = loop.start(); i <= loop.stop(); ++i ) {
			pose.set_phi( i, 50 );
			pose.set_psi( i, 100 );
		}

		CCDLoopClosureMover mover( loop, mm );
		mover.check_rama_scores( false );
		mover.max_per_move_torsion_delta_per_residue( 180., 180., 180. );
		mover.max_total_torsion_delta_per_residue( 180., 180., 180. );
		mover.max_cycles( 500 );
		mover.apply( pose );

		TS_ASSERT( mover.actual_cycles() < mover.max_cycles() );
		TS_ASSERT( mover.deviation() < mover.tolerance() );
		TS_ASSERT( mover.success() );
	}

	// test the copy constructor to make sure we are getting deep copies.
	void test_CCDLoopClosureMover_copy_ctor() {
		using namespace protocols::loops::loop_closure::ccd;

		CCDLoopClosureMover m1;

		// Set some values in m1 to non-defaults.
		m1.max_total_torsion_delta_per_residue( 100, 100, 100 );

		// initialize the Rama instance in m1
		RamaCheckBaseOP initial_rama = m1.rama();

		// copy m1
		CCDLoopClosureMover m2( m1 );

		// Make sure they start off the same
		TS_ASSERT( m1.max_total_torsion_delta_per_residue( helix ) == m2.max_total_torsion_delta_per_residue( helix ) );
		TS_ASSERT( m1.max_total_torsion_delta_per_residue( strand ) == m2.max_total_torsion_delta_per_residue( strand ) );
		TS_ASSERT( m1.max_total_torsion_delta_per_residue( coil ) == m2.max_total_torsion_delta_per_residue( coil ) );

		// change the values in m2
		m2.max_total_torsion_delta_per_residue( 0, 0, 0 );

		// Make sure they are different now
		TS_ASSERT( m1.max_total_torsion_delta_per_residue( helix ) != m2.max_total_torsion_delta_per_residue( helix ) );
		TS_ASSERT( m1.max_total_torsion_delta_per_residue( strand ) != m2.max_total_torsion_delta_per_residue( strand ) );
		TS_ASSERT( m1.max_total_torsion_delta_per_residue( coil ) != m2.max_total_torsion_delta_per_residue( coil ) );

		// Set m2 to use 2B Rama checks
		m2.use_rama_2B( true );
		TS_ASSERT( m2.rama()->name() == "RamaCheck2B" );
		TS_ASSERT( m1.rama()->name() == "RamaCheck1B" );

		// There are tons of checks that could be done here, but I'm satisifed for the time being

	}

	// test the MoverFactory to make sure we can create the mover and down cast appropriately
	void test_CCDLoopClosureMover_mover_factory() {
		using namespace protocols::moves;
		using namespace protocols::loops::loop_closure::ccd;

		std::string mover_name = "CCDLoopClosureMover";
		MoverFactory * mover_factory = MoverFactory::get_instance();
		MoverOP my_mover = mover_factory->newMover( mover_name );

		CCDLoopClosureMoverOP ccd_mover = utility::pointer::dynamic_pointer_cast< protocols::loops::loop_closure::ccd::CCDLoopClosureMover > ( my_mover );
		TS_ASSERT( ccd_mover ); // make sure we got back the right mover type

		// Create a Tag instance to test for RosettaScripts compatibility
		basic::datacache::DataMap dm;
		protocols::filters::Filters_map fm;
		Movers_map mm;

		using namespace utility::tag;
		using core::Real;

		// Set up tag so all options are not default
		TagOP tag( new Tag );
		tag->setName( mover_name );
		tag->setOption< Real >( "max_torsion_delta_per_move_H", 0 );
		tag->setOption< Real >( "max_torsion_delta_per_move_E", 0 );
		tag->setOption< Real >( "max_torsion_delta_per_move_L", 0 );

		tag->setOption< Real >( "max_torsion_delta_H", 0 );
		tag->setOption< Real >( "max_torsion_delta_E", 0 );
		tag->setOption< Real >( "max_torsion_delta_L", 0 );

		tag->setOption< Real >( "tolerance", 0 );

		tag->setOption< core::uint >( "max_cycles", 0 );
		tag->setOption< bool >( "check_rama_scores", false );
		tag->setOption< bool >( "rama_2b", true );

		// Only the Tag is used by the CCDLoopClosureMover's parse_my_tag, so the other instances that I'm passing through
		// can be nonsense
		MoverOP my_configured_mover = mover_factory->newMover( tag, dm, fm, mm, * pose_ );

		CCDLoopClosureMoverOP configured_ccd_mover = utility::pointer::dynamic_pointer_cast< protocols::loops::loop_closure::ccd::CCDLoopClosureMover > ( my_configured_mover );
		TS_ASSERT( configured_ccd_mover ); // make sure we got back the right mover type

		// Make sure the mover's configuration reflect what is in the Tag
		TS_ASSERT( configured_ccd_mover->max_per_move_torsion_delta_per_residue( helix ) == 0 );
		TS_ASSERT( configured_ccd_mover->max_per_move_torsion_delta_per_residue( strand ) == 0 );
		TS_ASSERT( configured_ccd_mover->max_per_move_torsion_delta_per_residue( coil ) == 0 );

		TS_ASSERT( configured_ccd_mover->max_total_torsion_delta_per_residue( helix ) == 0 );
		TS_ASSERT( configured_ccd_mover->max_total_torsion_delta_per_residue( strand ) == 0 );
		TS_ASSERT( configured_ccd_mover->max_total_torsion_delta_per_residue( coil ) == 0 );

		TS_ASSERT( configured_ccd_mover->tolerance() == 0 );

		TS_ASSERT( configured_ccd_mover->max_cycles() == 0 );
		TS_ASSERT( configured_ccd_mover->use_rama_2B() );
		TS_ASSERT( ! configured_ccd_mover->check_rama_scores() );
	}

	void test_CCDLoopClosureMover_exceptions() {
		Pose pose( *pose_ );
		CCDLoopClosureMover mover;

		// Test the two getters than throw key error exceptions.
		TR << "------------ A 'secstruct not valid' error message should follow -------------" << std::endl;
		try {
			set_throw_on_next_assertion_failure();
			mover.max_per_move_torsion_delta_per_residue( 'w' );
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			// ok -- let's break the error message into two lines.
			std::string msg = e.msg();
			std::vector< std::string > msg_lines = utility::split_by_newlines( msg );
			TS_ASSERT_EQUALS( msg_lines.size(), 3 );
			TS_ASSERT_EQUALS( msg_lines[1], "ERROR: CCDLoopClosureMover::max_per_move_delta_per_residue( char secstruct ): "
				"secstruct must be 'H', 'E', or 'L'. 'w' is not valid." );
		}
		TR << "------------ The previous error message was expected -------------" << std::endl;

		TR << "------------ A 'secstruct not valid' error message should follow -------------" << std::endl;
		try {
			set_throw_on_next_assertion_failure();
			mover.max_total_torsion_delta_per_residue( '!' );
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_Base const & e) {
			std::string msg = e.msg();
			std::vector< std::string > msg_lines = utility::split_by_newlines( msg );
			TS_ASSERT_EQUALS( msg_lines.size(), 3 );
			TS_ASSERT_EQUALS( msg_lines[1], "ERROR: CCDLoopClosureMover::max_total_delta_per_residue( char secstruct ): "
				"secstruct must be 'H', 'E', or 'L'. '!' is not valid." );
		}
		TR << "------------ The previous error message was expected -------------" << std::endl;

		// Test that a pose without cut points added will throw the user a message about that.
		// make sure we don't move omega
		MoveMapOP mm( new MoveMap );
		mm->set_bb( true );
		for ( core::uint i = 1; i < pose.total_residue(); ++i ) {
			mm->set( TorsionID( i, BB, 3 ), false );
		}

		Loop loop( 7, 14, 10 );

		set_single_loop_fold_tree( pose, loop );
		mover.movemap( mm );
		mover.loop( loop );
		TR << "------------ A 'Residue is not a cutpoint variant' error message should follow -------------" << std::endl;
		try {
			set_throw_on_next_assertion_failure();
			mover.apply( pose );
			TS_ASSERT( false );  // Exception was not thrown!
		} catch ( utility::excn::EXCN_BadInput const & e ) {
			TS_ASSERT_EQUALS( e.msg(), "CCDLoopClosureMover::get_anchors( core::conformation::Residue const & residue ): "
				"Residue is not a cutpoint variant! You must add cutpoint variants before applying this Mover." );
		}
		TR << "------------ The previous error message was expected -------------" << std::endl;

		// TODO: index_pair_in_range() is private, so I'm not sure how to test that one.... ~Labonte
	}
};

}  // anonymous namespace
