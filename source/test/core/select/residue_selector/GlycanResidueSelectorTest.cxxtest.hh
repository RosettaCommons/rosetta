// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/select/residue_selector/GlycanResidueSelectorTest.cxxtest.hh
/// @brief  Test glycan residue selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/core/select/residue_selector/utilities_for_testing.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/util.hh>

// Protocol Headers
#include <basic/Tracer.hh>

#include <utility/string_util.hh>

static basic::Tracer TR("GlycanResidueSelectorTest");

using namespace core::select;
using namespace core::select::residue_selector;
using core::Size;

class GlycanResidueSelectorTest : public CxxTest::TestSuite {



public:

	void setUp(){
		core_init_with_additional_options( "-include_sugars" );
		pose_from_file( pose_, "core/chemical/carbohydrates/gp120_2glycans_man5.pdb" , core::import_pose::PDB_file);



	}

	void test_residue_selectors(){

		///Vectors to test against
		utility::vector1< bool > correct_glycan_residues(pose_.size(), false);
		utility::vector1< bool > correct_single_branch_residues(pose_.size(), false);
		for ( core::Size i = 585; i <= 598; ++i ) {
			correct_glycan_residues[ i ] = true;
		}

		for ( core::Size i = 592; i <= 598; ++i ) {
			correct_single_branch_residues[ i ] = true;
		}

		//Test the selector

		GlycanResidueSelector selector = GlycanResidueSelector();


		utility::vector1< bool > glycan_residues = selector.apply(pose_);
		compare_bool_vector(correct_glycan_residues, glycan_residues);

		//Test get_residues_from_subset
		utility::vector1< Size > residues = get_residues_from_subset(glycan_residues);
		TS_ASSERT(residues.size() == 14);


		//Test getting residues of a branch
		utility::vector1< bool > root_residues(pose_.size(), false);
		root_residues[ 572 ] = true;
		selector.set_select_from_branch_residues(root_residues);
		glycan_residues = selector.apply(pose_);
		compare_bool_vector(correct_single_branch_residues, glycan_residues);



	}
	void tearDown(){

	}

private:

	core::pose::Pose pose_;







};



