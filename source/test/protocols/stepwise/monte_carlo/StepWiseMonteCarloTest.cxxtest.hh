// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/stepwise/monte_carlo//StepWiseMonteCarloTest.cxxtest.hh
/// @brief  test StepWiseMonteCarlo
/// @author Rhiju Das (rhiju@stanford.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocol Headers
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("StepWiseMonteCarloTest");


class StepWiseMonteCarloTest : public CxxTest::TestSuite {
	//Define Variables

public:

	protocols::stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo;
	core::pose::PoseOP pose;

	void setUp(){
		core_init_with_additional_options( "-fasta protocols/stepwise/modeler/align/seq_with_AAmismatch.fasta -s protocols/stepwise/modeler/align/scaff_subset_stepwise_input_seq_1.pdb" );

		using namespace core::chemical;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace protocols::stepwise::monte_carlo;
		using namespace protocols::stepwise::monte_carlo::options;
		using namespace protocols::stepwise::setup;

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		pose = initialize_pose_and_other_poses_from_command_line( rsd_set );

		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_res_level_energy.wts" ); // RNA/protein.
		stepwise_monte_carlo = StepWiseMonteCarloOP( new StepWiseMonteCarlo( scorefxn ) );
		StepWiseMonteCarloOptionsOP options( new StepWiseMonteCarloOptions );
		options->set_skip_preminimize( true );
		options->set_cycles( 0 ); // for bare bones run.
		stepwise_monte_carlo->set_options( options );
	}

	void tearDown(){
	}


	void test_skip_preminimize() {
		using namespace core::pose;

		PoseOP pose_save = pose->clone();
		stepwise_monte_carlo->apply( *pose );
		// better be a no op, since we set skip_preminimize and no cycles of monte carlo
		TS_ASSERT_EQUALS( pose->annotated_sequence(), pose_save->annotated_sequence() );
		for ( Size i = 1; i <= pose->total_residue(); ++i ) {
			for ( Size j = 1; j<= pose->residue_type( i ).natoms(); ++j ) {
				TS_ASSERT_LESS_THAN( ( pose->residue( i ).xyz( j ) - pose_save->residue( i ).xyz( j ) ).length(), 1.0e-6 );
			}
		}
	}

	void test_force_submotif_without_intervening_bulge(){
		using namespace protocols::stepwise::monte_carlo::options;

		TS_ASSERT( !stepwise_monte_carlo->master_mover()->stepwise_move_selector()->options()->force_submotif_without_intervening_bulge() );

		StepWiseMonteCarloOptionsOP options( stepwise_monte_carlo->options()->clone() );
		options->set_force_submotif_without_intervening_bulge( true );
		stepwise_monte_carlo->set_options( options );
		TS_ASSERT( stepwise_monte_carlo->master_mover()->stepwise_move_selector()->options()->force_submotif_without_intervening_bulge() );
	}



};



