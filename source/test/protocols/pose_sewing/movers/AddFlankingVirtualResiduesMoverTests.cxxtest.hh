// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pose_sewing/movers/AddFlankingVirtualResiduesMoverTests.cxxtest.hh
/// @brief  tests the AddFlankingVirtualResiduesMover
/// @author Frank Teets (frankdt@email.unc.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/pose_sewing/movers/AddFlankingVirtualResiduesMover.hh>
#include <core/select/residue_selector/SecondaryStructureSelector.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("AddFlankingVirtualResiduesMoverTests");


class AddFlankingVirtualResiduesMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}
	void test_AddFlankingVirtualResiduesMover(){
		core::select::residue_selector::SecondaryStructureSelectorOP in_selector( new core::select::residue_selector::SecondaryStructureSelector );
		in_selector->set_selected_ss("L");
		in_selector->set_use_dssp(true);
		protocols::pose_sewing::movers::AddFlankingVirtualResiduesMoverOP mover(new protocols::pose_sewing::movers::AddFlankingVirtualResiduesMover);
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::Size orig_size = trpcage.size();
		mover->set_N_term_length(10);
		mover->set_C_term_length(10);
		mover->set_chain_to_modify(1);
		mover->set_vital_selector(in_selector);
		mover->add_flanking_virtual_residues(trpcage);
		TS_ASSERT(trpcage.size() == orig_size+20);

	}

	void tearDown() {

	}






};
