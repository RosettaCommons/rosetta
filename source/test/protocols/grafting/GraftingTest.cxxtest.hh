// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/grafting/GraftingTest.cxxtest.hh
/// @brief  tests for container GraftMovers and utility functions
/// @author Jared Adolf-Bryfogle, Brian Weitzner

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/grafting/simple_movers/InsertPoseIntoPoseMover.hh>
#include <protocols/grafting/simple_movers/ReplaceRegionMover.hh>
#include <protocols/grafting/simple_movers/KeepRegionMover.hh>

#include <protocols/grafting/util.hh>

// Core Headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

// Protocol Headers
#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.grafting.GraftTest");
class GraftingTest : public CxxTest::TestSuite {
	core::pose::Pose scaffold_pose; //Full PDB
	core::pose::Pose piece; //CDR to graft.

	core::Size start;
	core::Size end;
	core::Size flex;
	core::Size nter_overhang;
	core::Size cter_overhang;
	core::Size starting_residues;
	core::Size insert_size;
	core::Size L1_2j88_size;


public:

	void setUp(){

		core_init();
		core::import_pose::pose_from_file(scaffold_pose, "protocols/antibody/2j88.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file(piece, "protocols/grafting/2aabL_L1.pdb", core::import_pose::PDB_file);

		starting_residues = scaffold_pose.size();
		nter_overhang=3;
		cter_overhang=3;
		flex=2;
		start = scaffold_pose.pdb_info()->pdb2pose('L', 24) - 1;
		end = scaffold_pose.pdb_info()->pdb2pose('L', 42) + 1;
		insert_size = piece.size()-nter_overhang-cter_overhang;
		L1_2j88_size = 11;

		TR <<"Setup"<<std::endl;
	}

	void tearDown(){
		scaffold_pose.clear();
		piece.clear();
	}

	void test_utility_functions(){
		using namespace core::kinematics;
		using namespace protocols::grafting::simple_movers;
		TR<<"Return region"<<std::endl;
		core::pose::Pose new_region = protocols::grafting::return_region(scaffold_pose, scaffold_pose.pdb_info()->pdb2pose('L', 24), scaffold_pose.pdb_info()->pdb2pose('L', 42));
		TS_ASSERT_EQUALS(L1_2j88_size, new_region.size());

		TR<<"Replace Region"<<std::endl;
		TR << new_region << std::endl;
		core::pose::Pose replaced_pose = protocols::grafting::replace_region(new_region, scaffold_pose, 1, 24, new_region.size(), true);
		TS_ASSERT_EQUALS(starting_residues, replaced_pose.size());

		TR <<"Replace Region Mover" << std::endl;
		TR << new_region << std::endl;
		core::pose::Pose scaffold_copy = core::pose::Pose(scaffold_pose);
		ReplaceRegionMover replacer = ReplaceRegionMover(new_region, 1, 24, new_region.size());
		replacer.apply(scaffold_copy);
		TS_ASSERT_EQUALS(starting_residues, scaffold_pose.size());

		TR<<"Delete Region"<<std::endl;
		scaffold_copy = core::pose::Pose(scaffold_pose);
		protocols::grafting::delete_region(scaffold_copy, start+1, end-1);
		core::Size deleted_residues  = starting_residues-scaffold_copy.size();
		TS_ASSERT_EQUALS(deleted_residues, L1_2j88_size);

		TR<<"Insert Region"<<std::endl;
		core::pose::Pose new_pose = protocols::grafting::insert_pose_into_pose(scaffold_copy, new_region, start, start+1, true);
		TS_ASSERT_EQUALS(new_pose.size(), starting_residues);

		TR<<"Delete Region Mover" << std::endl;
		scaffold_copy = core::pose::Pose(scaffold_pose);
		DeleteRegionMover deleter = DeleteRegionMover(start+1, end-1);
		deleter.apply(scaffold_copy);
		deleted_residues  = starting_residues-scaffold_copy.size();
		TS_ASSERT_EQUALS(deleted_residues, L1_2j88_size);

		TR<<"Insert Region Mover" << std::endl;
		InsertPoseIntoPoseMover inserter = InsertPoseIntoPoseMover(new_region, start, start+1);
		inserter.apply(scaffold_copy);
		TS_ASSERT_EQUALS(new_pose.size(), starting_residues);

		TR << "Keep Region Mover" << std::endl;
		scaffold_copy = core::pose::Pose(scaffold_pose);
		KeepRegionMover keeper = KeepRegionMover(start+1, end-1);
		keeper.apply(scaffold_copy);
		TS_ASSERT_EQUALS(scaffold_copy.size(), L1_2j88_size);

		TR << "Combine MoveMaps" << std::endl;

		MoveMapOP scaffold_mm( new MoveMap() );
		MoveMapOP insert_mm( new MoveMap() );

		scaffold_mm->set_bb(22, true);
		scaffold_mm->set_bb(23, true);
		scaffold_mm->set_bb(35, true);
		scaffold_mm->set_bb(36, true);

		insert_mm->set_bb(nter_overhang+1, true);
		insert_mm->set_bb(nter_overhang+2, true);
		insert_mm->set_bb(piece.size()-cter_overhang, true);
		insert_mm->set_bb(piece.size()-cter_overhang - 1, true);

		MoveMapOP combined_mm = protocols::grafting::combine_movemaps_post_insertion(
			scaffold_mm, insert_mm, start, end, piece.size() - nter_overhang - cter_overhang, cter_overhang);

		for ( core::Size i = 1; i <= 21; ++i ) {
			TS_ASSERT(combined_mm->get_bb(i) == false);
		}

		//Either two residues at either side of Nter graft connection
		TS_ASSERT(combined_mm->get_bb(22) == true);
		TS_ASSERT(combined_mm->get_bb(23) == true);
		TS_ASSERT(combined_mm->get_bb(24) == true);
		TS_ASSERT(combined_mm->get_bb(25) == true);

		//Either two residues at either side of Cter graft connection
		TS_ASSERT(combined_mm->get_bb(37) == true);
		TS_ASSERT(combined_mm->get_bb(38) == true);
		TS_ASSERT(combined_mm->get_bb(39) == true);
		TS_ASSERT(combined_mm->get_bb(40) == true);

		for ( core::Size i = 41; i <=starting_residues + 4; ++i ) {
			TS_ASSERT(combined_mm->get_bb(i) == false);
		}
		TR << "Combining Movemaps Complete." << std::endl;
	}

	void test_graft_classes(){
		using namespace protocols::grafting;


		//////////////CCD Ends Graft ///////////////////////////////////
		TR << "Testing CCDEndsGraftMover Mover" << std::endl;
		TR << "start: " << start << std::endl;
		TR << "end: " << end << std::endl;

		CCDEndsGraftMoverOP cdr_grafter( new protocols::grafting::CCDEndsGraftMover(start, end, piece, nter_overhang, cter_overhang) );
		core::pose::Pose scaffold_copy = core::pose::Pose(scaffold_pose);

		cdr_grafter->set_scaffold_flexibility(flex, flex);
		cdr_grafter->set_skip_sampling(false);

		cdr_grafter->final_repack(true);
		cdr_grafter->stop_at_closure(true);

		cdr_grafter->set_insert_flexibility(0, 0);
		cdr_grafter->set_cycles(5);
		cdr_grafter->apply(scaffold_copy);

		TS_ASSERT_EQUALS(scaffold_copy.size(), starting_residues - L1_2j88_size + insert_size);
		TS_ASSERT_EQUALS(cdr_grafter->start(), start);
		TS_ASSERT_EQUALS(cdr_grafter->original_end(), end);
		TS_ASSERT_EQUALS(cdr_grafter->end(), 39);
		TS_ASSERT_EQUALS(cdr_grafter->insertion_length(), insert_size);


		////////////Anchored Graft /////////////////////////////////////
		TR<< "Testing AnchoredGraft Mover" << std::endl;
		AnchoredGraftMoverOP anchored_grafter( new protocols::grafting::AnchoredGraftMover(start, end, true) );
		scaffold_copy = core::pose::Pose(scaffold_pose);

		anchored_grafter->set_piece(piece, nter_overhang, cter_overhang);
		anchored_grafter->set_scaffold_flexibility(flex, flex);
		anchored_grafter->set_insert_flexibility(0, 0);
		anchored_grafter->final_repack(false);
		anchored_grafter->stop_at_closure(false);
		anchored_grafter->set_cycles(2);
		anchored_grafter->copy_pdbinfo(true);
		anchored_grafter->apply(scaffold_copy);

		TS_ASSERT_EQUALS(scaffold_copy.size(), starting_residues - L1_2j88_size + insert_size);
		TS_ASSERT_EQUALS(anchored_grafter->start(), start);
		TS_ASSERT_EQUALS(anchored_grafter->original_end(), end);
		TS_ASSERT_EQUALS(anchored_grafter->end(), 39); //Same length
		TS_ASSERT_EQUALS(anchored_grafter->insertion_length(), insert_size);

		//Check PDBInfo Copy - first residue of L1 from insert start.
		// If copy failed - we cannot access.
		TS_ASSERT_EQUALS(scaffold_copy.pdb_info()->pdb2pose('L', 23), start);
		TS_ASSERT_EQUALS(scaffold_copy.pdb_info()->pdb2pose('L', 24), start + 1);

	}


};


