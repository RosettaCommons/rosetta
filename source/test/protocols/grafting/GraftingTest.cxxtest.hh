// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

// Protocol Headers
#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.grafting.GraftTest");
class GraftingTest : public CxxTest::TestSuite {
	core::pose::Pose scaffold_pose; //Full PDB
	core::pose::Pose framework_pose; //PDB Missing a cdr.
	core::pose::Pose piece; //CDR to graft.
	
	core::Size start;
	core::Size end;
	core::Size flex;
	core::Size nter_overhang;
	core::Size cter_overhang;
	core::Size starting_residues;
	core::Size insert_size;

    
public:
	
	void setUp(){
		
		core_init();
		core::import_pose::pose_from_pdb(scaffold_pose, "protocols/antibody/2j88.pdb");
		core::import_pose::pose_from_pdb(framework_pose, "protocols/grafting/2j88_NoL1.pdb");
		core::import_pose::pose_from_pdb(piece, "protocols/grafting/2j88_L1_overhang3_rotated.pdb");
		
		starting_residues = scaffold_pose.total_residue();
		nter_overhang=3;
		cter_overhang=3;
		flex=2;
		start = 23;
		end = 35;
		insert_size = piece.total_residue()-nter_overhang-cter_overhang;
		TR <<"Setup"<<std::endl;
	}
	
	void tearDown(){
		scaffold_pose.clear();
		framework_pose.clear();
		piece.clear();
	}
    
	void test_utility_functions(){
		using namespace core::kinematics;
		using namespace protocols::grafting::simple_movers;
		TR<<"Return region"<<std::endl;
		core::pose::Pose new_region = protocols::grafting::return_region(piece, 1+nter_overhang, piece.total_residue()-cter_overhang);
		TS_ASSERT_EQUALS(insert_size, new_region.total_residue());
		
		TR<<"Replace Region"<<std::endl;
		TR << new_region << std::endl;
		core::pose::Pose replaced_pose = protocols::grafting::replace_region(new_region, scaffold_pose, 1, 24, new_region.total_residue(), true);
		TS_ASSERT_EQUALS(starting_residues, replaced_pose.total_residue());

		TR <<"Replace Region Mover" << std::endl;
		TR << new_region << std::endl;
		core::pose::Pose scaffold_copy = core::pose::Pose(scaffold_pose);
		ReplaceRegionMover replacer = ReplaceRegionMover(new_region, 1, 24, new_region.total_residue());
		replacer.apply(scaffold_copy);
		TS_ASSERT_EQUALS(starting_residues, scaffold_pose.total_residue());
		
		TR<<"Delete Region"<<std::endl;
		scaffold_copy = core::pose::Pose(scaffold_pose);
		protocols::grafting::delete_region(scaffold_copy, start+1, end-1);
		core::Size deleted_residues  = starting_residues-scaffold_copy.total_residue();
		TS_ASSERT_EQUALS(deleted_residues, insert_size);
	
		TR<<"Insert Region"<<std::endl;
		core::pose::Pose new_pose = protocols::grafting::insert_pose_into_pose(scaffold_copy, new_region, start, start+1, true);
		TS_ASSERT_EQUALS(new_pose.total_residue(), starting_residues);
		
		TR<<"Delete Region Mover" << std::endl;
		scaffold_copy = core::pose::Pose(scaffold_pose);
		DeleteRegionMover deleter = DeleteRegionMover(start+1, end-1);
		deleter.apply(scaffold_copy);
		deleted_residues  = starting_residues-scaffold_copy.total_residue();
		TS_ASSERT_EQUALS(deleted_residues, insert_size);
		
		TR<<"Insert Region Mover" << std::endl;
		InsertPoseIntoPoseMover inserter = InsertPoseIntoPoseMover(new_region, start, start+1);
		inserter.apply(scaffold_copy);
		TS_ASSERT_EQUALS(new_pose.total_residue(), starting_residues);
		
		TR << "Keep Region Mover" << std::endl;
		scaffold_copy = core::pose::Pose(scaffold_pose);
		KeepRegionMover keeper = KeepRegionMover(start+1, end-1);
		keeper.apply(scaffold_copy);
		TS_ASSERT_EQUALS(scaffold_copy.total_residue(), insert_size);
		
		TR << "Combine MoveMaps" << std::endl;
		
		MoveMapOP scaffold_mm = new MoveMap();
		MoveMapOP insert_mm = new MoveMap();
		
		scaffold_mm->set_bb(23, true);
		scaffold_mm->set_bb(22, true);
		scaffold_mm->set_bb(35, true);
		scaffold_mm->set_bb(36, true);
		
		insert_mm->set_bb(nter_overhang+1, true);
		insert_mm->set_bb(nter_overhang+2, true);
		insert_mm->set_bb(piece.total_residue()-cter_overhang, true);
		insert_mm->set_bb(piece.total_residue()-cter_overhang - 1, true);
		
		MoveMapOP combined_mm = protocols::grafting::combine_movemaps_post_insertion(
			scaffold_mm, insert_mm, 23, 35, piece.total_residue() - nter_overhang - cter_overhang, cter_overhang);
		
		TS_ASSERT(combined_mm->get_bb(21) == false);
		TS_ASSERT(combined_mm->get_bb(23) == true);
		TS_ASSERT(combined_mm->get_bb(22) == true);
		TS_ASSERT(combined_mm->get_bb(35) == true);
		TS_ASSERT(combined_mm->get_bb(36) == true);
		TS_ASSERT(combined_mm->get_bb(37) == false);
		TS_ASSERT(combined_mm->get_bb(24) == true);
		TS_ASSERT(combined_mm->get_bb(25) == true);
		TS_ASSERT(combined_mm->get_bb(34) == true);
		TS_ASSERT(combined_mm->get_bb(33) == true);
		
	}
    
	void test_graft_classes(){
		using namespace protocols::grafting;
		AnchoredGraftMoverOP anchored_grafter = new protocols::grafting::AnchoredGraftMover(start, end, true);
		CCDEndsGraftMoverOP cdr_grafter = new protocols::grafting::CCDEndsGraftMover(start, end, piece, nter_overhang, cter_overhang);
		
		////////////Anchored Graft /////////////////////////////////////
		core::pose::Pose scaffold_copy = core::pose::Pose(scaffold_pose);
		
		
		//core::Real rms = core::scoring::CA_rmsd(scaffold_copy, scaffold_pose);
		//core::Real rms_assert = 0.0;
		//TS_ASSERT_EQUALS(rms, rms_assert);
		
		//////////////CCD Ends Graft ///////////////////////////////////
		scaffold_copy = core::pose::Pose(scaffold_pose);
		TR<<"Testing CCDEndsGraftMover Mover"<<std::endl;
		cdr_grafter->set_scaffold_flexibility(flex, flex);
		cdr_grafter->set_skip_sampling(false);
		
		cdr_grafter->final_repack(true);
		cdr_grafter->stop_at_closure(true);
		
		cdr_grafter->set_insert_flexibility(0, 0);
		cdr_grafter->set_cycles(5);
		cdr_grafter->apply(scaffold_copy);
		//scaffold_copy.dump_pdb("/home/jadolfbr/Documents/modeling/rosetta/Rosetta/main/source/test.pdb");
		TS_ASSERT_EQUALS(scaffold_copy.total_residue(), starting_residues);
		//core::Real rms = core::scoring::CA_rmsd(scaffold_copy, scaffold_pose);
		//core::Real rms_assert = 0.0;
		
		
		TR<< "Testing AnchoredGraft Mover" << std::endl;
		
		anchored_grafter->set_piece(piece, nter_overhang, cter_overhang);
		
		anchored_grafter->set_scaffold_flexibility(flex, flex);
		anchored_grafter->set_insert_flexibility(0, 0);
		anchored_grafter->final_repack(false);
		anchored_grafter->stop_at_closure(false);
		anchored_grafter->set_cycles(2);
		anchored_grafter->apply(scaffold_copy);
		//scaffold_copy.dump_pdb("/home/jadolfbr/Documents/modeling/rosetta/Rosetta/main/source/test2.pdb");
		
		TS_ASSERT_EQUALS(scaffold_pose.total_residue(), starting_residues);
		
	}
    
};



