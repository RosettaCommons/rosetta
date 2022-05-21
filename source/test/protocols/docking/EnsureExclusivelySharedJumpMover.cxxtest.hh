// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/docking/EnsureExclusivelySharedJumpMover.cxxtest.hh
/// @brief The asserts themselves are inside EnsureExclusivelySharedJumpMover::apply. These tests just set up different problems to pass in to those asserts.
/// @author Jack Maguire

// Test headers
#include <test/UTracer.hh>

// Unit headers
#include <protocols/docking/EnsureExclusivelySharedJumpMover.hh>
#include <core/select/jump_selector/ExclusivelySharedJumpSelector.hh>

// Package headers
#include <protocols/docking/metrics.hh>
#include <protocols/docking/util.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh> // AUTO IWYU For pose_from_file, PDB_file
#include <core/init_util.hh> // AUTO IWYU For core_init_with_additional_options
#include <core/kinematics/MoveMap.hh> // AUTO IWYU For MoveMap
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <basic/Tracer.hh>

using namespace protocols;
using namespace protocols::docking;
using namespace core::select::residue_selector;
using namespace core::select::jump_selector;
using namespace core::kinematics;

static basic::Tracer TR( "tests.protocols.docking.EnsureExclusivelySharedJumpMoverTests" );

class EnsureExclusivelySharedJumpMoverTests : public CxxTest::TestSuite {
public:
	void setUp() {
		core_init();
	}

	void test_simple(){
		core::pose::Pose pose;
		std::string const seq(20, 'A');
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "1-20" ) );

		//This should throw because everything is selected
		TS_ASSERT_THROWS_ANYTHING( mover.apply( pose ) );
	}

	void test_simple_two_chain1(){
		core::pose::Pose pose;
		std::string const seq = "AAA/GGG";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "1-3" ) );

		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );
	}

	void test_simple_two_chain2(){
		core::pose::Pose pose;
		std::string const seq = "AAA/GGG";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "4-6" ) );

		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );
	}

	void test_simple_two_chain3(){
		core::pose::Pose pose;
		std::string const seq = "AAA/GGG";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "2-3" ) );

		//This should throw because residues 1 and 2 cannot be split by a jump
		TS_ASSERT_THROWS_ANYTHING( mover.apply( pose ) );
	}

	void test_simple_two_chain4(){
		core::pose::Pose pose;
		std::string const seq = "AAA/GGG";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "4,6" ) );

		//This should throw because residues 4 and 5 cannot be split by a jump
		TS_ASSERT_THROWS_ANYTHING( mover.apply( pose ) );
	}


	void test_custom_foldtree1(){
		core::pose::Pose pose;
		std::string const seq = "AAA/GGG/PPP";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "4-6" ) );

		//This is fine
		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );

		{
			FoldTree ft = pose.fold_tree();
			ft.reorder( 4 );
			pose.fold_tree( ft );
		}

		//This is fine
		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );

		{
			FoldTree ft = pose.fold_tree();
			ft.reorder( 7 );
			pose.fold_tree( ft );
		}

		//This is fine
		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );
	}

	void test_custom_foldtree2(){
		core::pose::Pose pose;
		std::string const seq = "AAA/GGG/PPP";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "4-9" ) );

		//This is fine
		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );

		{
			FoldTree ft = pose.fold_tree();
			ft.reorder( 4 );
			pose.fold_tree( ft );
		}

		//This is fine
		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );

		{
			FoldTree ft = pose.fold_tree();
			ft.reorder( 7 );
			pose.fold_tree( ft );
		}

		//This is fine
		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );
	}

	void test_custom_foldtree3(){
		core::pose::Pose pose;
		std::string const seq = "AAAGGGPPP";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		auto const sele = utility::pointer::make_shared< ResidueIndexSelector >( "4-6" );
		EnsureExclusivelySharedJumpMover mover( sele );

		FoldTree ft( pose.size() );
		ft.clear();
		ft.add_edge( 1, 3, -1 );
		ft.add_edge( 1, 4, 1 );
		ft.add_edge( 4, 6, -1 );
		ft.add_edge( 1, 7, 2 );
		ft.add_edge( 7, 9, -1 );
		pose.fold_tree( ft );

		//This is fine
		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );

		ExclusivelySharedJumpSelector const js( sele );
		auto const jumps = js.selection_jumps( pose );
		TS_ASSERT_EQUALS( jumps.size(), 1 );
		TS_ASSERT_EQUALS( jumps[1], 1 );
	}

	void test_custom_foldtree4(){
		core::pose::Pose pose;
		std::string const seq = "AAAGGGPPP";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		auto sele = utility::pointer::make_shared< ResidueIndexSelector >( "4-9" );
		EnsureExclusivelySharedJumpMover mover( sele );

		FoldTree ft( pose.size() );
		ft.clear();
		ft.add_edge( 1, 3, -1 );
		ft.add_edge( 1, 4, 1 );
		ft.add_edge( 4, 6, -1 );
		ft.add_edge( 1, 7, 2 );
		ft.add_edge( 7, 9, -1 );
		pose.fold_tree( ft );

		//This is fine
		TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );

		ExclusivelySharedJumpSelector const js( sele );
		auto const jumps = js.selection_jumps( pose );
		for ( auto const & e : pose.fold_tree() ) {
			TR << "e: " << e << std::endl;
		}
		for ( core::Size const j : jumps ) {
			TR << "j: " << j << std::endl;
		}
		TS_ASSERT_EQUALS( jumps.size(), 1 );
		TS_ASSERT_EQUALS( jumps[1], 1 );
	}

	void test_custom_foldtree5(){
		core::pose::Pose pose;
		std::string const seq = "AAAGGGPPP";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");
		EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "4,8,9" ) );

		FoldTree ft( pose.size() );
		ft.clear();
		ft.add_edge( 1, 3, -1 );
		ft.add_edge( 1, 4, 1 );
		ft.add_edge( 4, 6, -1 );
		ft.add_edge( 1, 7, 2 );
		ft.add_edge( 7, 9, -1 );
		pose.fold_tree( ft );

		//This is not fine
		TS_ASSERT_THROWS_ANYTHING( mover.apply( pose ) );
	}

	void test_simple_six_chains(){
		core::pose::Pose pose;
		std::string const seq = "AAA/AAA/AAA/AAA/AAA/AAA";
		core::pose::make_pose_from_sequence(pose, seq, "fa_standard");

		{
			TR << "1-3" << std::endl;
			EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "1-3" ) );
			TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );
		}

		{
			TR << "1-18" << std::endl;
			EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "1-18" ) );
			//Can't have everything selected
			TS_ASSERT_THROWS_ANYTHING( mover.apply( pose ) );
		}

		{
			TR << "1-6" << std::endl;
			EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "1-6" ) );
			TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );
		}

		{
			TR << "1-7" << std::endl;
			EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "1-7" ) );
			//Can't split chains
			TS_ASSERT_THROWS_ANYTHING( mover.apply( pose ) );
		}

		{
			TR << "1-15" << std::endl;
			EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "1-15" ) );
			TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );
		}

		{
			TR << "1-3,16-18" << std::endl;
			EnsureExclusivelySharedJumpMover mover( utility::pointer::make_shared< ResidueIndexSelector >( "1-3,16-18" ) );
			TS_ASSERT_THROWS_NOTHING( mover.apply( pose ) );
		}

	}

};
