// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.cxxtest.hh.
/// @brief  Unit tests for the simple_cycpep_predict application.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers:
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pdb1ubq.hh>

// GeneralizedKIC headers:
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI_JobResultsSummary.hh>
#include <protocols/cyclic_peptide_predict/util.hh>

// Other Rosetta libraries:
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/MutateResidue.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.cyclic_peptide_predict.simple_cycpep_predict_design.cxxtest.hh");


// --------------- Test Class --------------- //

class SimpleCycpepPredictApplication_design_Tests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP testpose_;
	core::scoring::ScoreFunctionOP scorefxn_;

public:

	void setUp() {
		core_init_with_additional_options( "-cyclic_peptide:sequence_file protocols/cyclic_peptide_predict/seq_gly.txt -cyclic_peptide:allowed_residues_by_position protocols/cyclic_peptide_predict/allowed_residues.txt -cyclic_peptide:design_peptide -nstruct 1 -symmetric_gly_tables true -cyclic_peptide:genkic_closure_attempts 1000 -cyclic_peptide:genkic_min_solution_count 1 -cyclic_peptide:min_genkic_hbonds 1 -cyclic_peptide:min_final_hbonds 1 -cyclic_peptide:fast_relax_rounds 1 -cyclic_peptide:rama_cutoff 3.0 -in:file:native protocols/cyclic_peptide_predict/native.pdb" );
	}

	void tearDown() {
	}

	/// @brief This is the simplest test in the world: run the protocol and see if it crashes.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_simple_cycpep_predict() {
		TR << "Running SimpleCycpepPredictApplication_design_Tests::test_simple_cycpep_predict()." << std::endl;
		TR << "This simply runs the SimpleCycpepPredictApplication in design mode and checks for crashes." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;
		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication the_app(true);
		the_app.run();
	}

	/// @brief This is the simplest test in the world: run the protocol and see if it crashes.  This
	/// version uses a setup file.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_simple_cycpep_predict_withfile() {
		TR << "Running SimpleCycpepPredictApplication_design_Tests::test_simple_cycpep_predict_withfile()." << std::endl;
		TR << "This simply runs the SimpleCycpepPredictApplication in design mode and checks for crashes." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;
		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication the_app;
		the_app.run();
	}

}; //class GeneralizedKIC_Tests
