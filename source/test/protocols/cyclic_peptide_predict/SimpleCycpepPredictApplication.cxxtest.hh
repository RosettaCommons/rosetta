// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.cxxtest.hh.
/// @brief  Unit tests for the simple_cycpep_predict application.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers:
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pdb1ubq.hh>

// GeneralizedKIC headers:
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>

// Other Rosetta libraries:
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/MutateResidue.hh>


// --------------- Test Class --------------- //

class SimpleCycpepPredictApplicationTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
	core::scoring::ScoreFunctionOP scorefxn_;

public:

	void setUp() {
		core_init_with_additional_options( "-cyclic_peptide:sequence_file protocols/cyclic_peptide_predict/seq.txt -symmetric_gly_tables true -cyclic_peptide:genkic_closure_attempts 1000 -cyclic_peptide:genkic_min_solution_count 1 -cyclic_peptide:min_genkic_hbonds 2 -cyclic_peptide:min_final_hbonds 2 -cyclic_peptide:fast_relax_rounds 1 -cyclic_peptide:rama_cutoff 3.0 -in:file:native protocols/cyclic_peptide_predict/native.pdb" );
	}

	void tearDown() {
	}

	/// @brief This is shte simplest test in the world: run the protocol and see if it crashes.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_simple_cycpep_predict() {
		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication the_app;
		the_app.run();
	}


}; //class GeneralizedKIC_Tests
