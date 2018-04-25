// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/select/residue_selector/ResidueInSequenceMotifSelectorTest.cxxtest.hh
/// @brief  Test sequence motif and GlycanSequons selectors using regular expressions
/// @author Sebastian Rämsch (raemisch@scripps.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/core/select/residue_selector/utilities_for_testing.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/select/residue_selector/ResidueInSequenceMotifSelector.hh>
#include <core/select/residue_selector/GlycanSequonsSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

// Protocol Headers
#include <basic/Tracer.hh>

#include <utility/string_util.hh>

static basic::Tracer TR("ResidueInSequenceMotifSelectorTest");

using namespace core::select;
using namespace core::select::residue_selector;
using core::Size;

class ResidueInSequenceMotifSelectorTest : public CxxTest::TestSuite {

public:

	void setUp(){
		core_init_with_additional_options( "-include_sugars" );
		pose_from_file( pose_, "core/chemical/carbohydrates/glycosylated_peptide.pdb" , core::import_pose::PDB_file);
	}

	void test_residue_selectors(){

		// this cannot be run on systems with "older" versions of the c++11 std library
		// those without regex support will compile and link but not run.
		if ( ! regex_usable() ) { return; }

		///Vector to test against
		utility::vector1< bool > correct_residues(pose_.size(), false);
		correct_residues[ 3 ] = true;

		// Sequence of pose: "HLNSS"

		// Test the selector
		// Find NxS or NxT
		ResidueInSequenceMotifSelector selector = ResidueInSequenceMotifSelector( "N.[ST]",1 );
		utility::vector1< bool > selected_residues = selector.apply(pose_);
		compare_bool_vector(correct_residues, selected_residues);

		GlycanSequonsSelector sequon_selector = GlycanSequonsSelector();
		selected_residues = sequon_selector.apply(pose_);
		compare_bool_vector(correct_residues, selected_residues);

		// test for multiple occurances
		// here: find either N or S
		// Sequence of pose: "HLNSS"
		correct_residues[ 4 ] = true;
		correct_residues[ 5 ] = true;

		selector = ResidueInSequenceMotifSelector( "[NS]",1 );
		selected_residues = selector.apply(pose_);
		compare_bool_vector(correct_residues, selected_residues);

	}
	void tearDown(){

	}

private:

	core::pose::Pose pose_;
};
