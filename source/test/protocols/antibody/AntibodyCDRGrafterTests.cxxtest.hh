// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/AntibodyCDRGrafterTests
/// @brief  tests for container GraftMovers and utility functions
/// @author Jared Adolf-Bryfogle, Brian Weitzner

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/grafting/simple_movers/InsertPoseIntoPoseMover.hh>
#include <protocols/grafting/simple_movers/ReplaceRegionMover.hh>
#include <protocols/grafting/simple_movers/KeepRegionMover.hh>

#include <protocols/antibody/AntibodyCDRGrafter.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/grafting/util.hh>

// Core Headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

// Protocol Headers
#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.antibody.AntibodyCDRGrafterTests");
class AntibodyCDRGrafterTests : public CxxTest::TestSuite {
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

	void test_graft_classes(){
		using namespace protocols::grafting;
		using namespace protocols::antibody;


		//////////////CCD Ends Graft ///////////////////////////////////
		TR << "Testing AntibodyCDRGrafter Mover" << std::endl;
		TR << "start: " << start << std::endl;
		TR << "end: " << end << std::endl;
		TR <<"nter overhang: " << nter_overhang << std::endl;
		TR <<"cter overhang: " << cter_overhang << std::endl;
		TR << "Flex: " << flex << std::endl;
		TR << "Insert size" << insert_size << std::endl;

		AntibodyInfoOP ab_info( new AntibodyInfo(scaffold_pose, AHO_Scheme, North));
		TR << "Created Ab info for scaffold" << std::endl;

		AntibodyCDRGrafterOP grafter( new AntibodyCDRGrafter(ab_info));

		TR << "Set grafter up" << std::endl;
		grafter->set_cdr_only(l1);
		grafter->set_donor_structure(piece);

		//Test some functions
		grafter->set_stop_after_closure(true);
		grafter->set_use_secondary_graft_mover_if_needed(false);
		grafter->set_idealize_insert(true);

		core::pose::Pose scaffold_copy = core::pose::Pose(scaffold_pose);

		TR << "Ready to roll!" << std::endl;
		grafter->apply(scaffold_copy);

		TR << "Graft complete" << std::endl;

		TS_ASSERT_EQUALS(scaffold_copy.size(), starting_residues - L1_2j88_size + insert_size);
		TS_ASSERT_EQUALS(grafter->get_primary_graft_mover().start(), start);
		TS_ASSERT_EQUALS(grafter->get_primary_graft_mover().original_end(), end);
		TS_ASSERT_EQUALS(grafter->get_primary_graft_mover().end(), 39);
		TS_ASSERT_EQUALS(grafter->get_primary_graft_mover().insertion_length(), insert_size);


		//Check PDBInfo Copy - first residue of L1 from insert start.
		// If copy failed - we cannot access.
		TS_ASSERT_EQUALS(scaffold_copy.pdb_info()->pdb2pose('L', 23), start);
		TS_ASSERT_EQUALS(scaffold_copy.pdb_info()->pdb2pose('L', 24), start + 1);

	}


};


