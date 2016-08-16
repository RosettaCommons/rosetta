// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/stepwise/modeler//StepWiseMinimizerTest.cxxtest.hh
/// @brief  test StepWiseMinimizer
/// @author Rhiju Das (rhiju@stanford.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/FoldTree.hh>

// Protocol Headers
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/stepwise/modeler/StepWiseMinimizer.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/packer/StepWisePacker.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/working_parameters/util.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.hh>
#include <protocols/stepwise/monte_carlo/util.hh>

#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("StepWiseMinimizerTest");


class StepWiseMinimizerTest : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		// same setup as stepwise/modeler/align/StepWisePoseAlignerTest.cxxtest.hh
		core_init_with_additional_options( "-s protocols/stepwise/modeler/align/1_scaff_subset_stepwise_input_seq_1.pdb protocols/stepwise/modeler/align/2_scaff_subset_stepwise_input_seq_1.pdb  -fasta protocols/stepwise/modeler/align/seq_with_AAmismatch.fasta -align_pdb protocols/stepwise/modeler/align/scaff_subset_stepwise_input_seq_1.pdb" );
	}

	void tearDown(){
	}

	void test_num_pose_minimize(){

		using namespace core::chemical;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		using namespace core::id;
		using namespace protocols::stepwise::setup;
		using namespace protocols::stepwise::monte_carlo;
		using namespace protocols::stepwise::monte_carlo::mover;
		using namespace protocols::stepwise::modeler;
		using namespace protocols::stepwise::modeler::options;
		using namespace protocols::stepwise::modeler::working_parameters;
		using namespace protocols::stepwise::modeler::packer;
		using namespace utility;
		using namespace utility::tools;

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		PoseOP pose_op = initialize_pose_and_other_poses_from_command_line( rsd_set );
		pose::Pose & pose = *pose_op;

		AddMover add_mover;
		add_mover.set_start_added_residue_in_aform( true );
		add_mover.set_presample_added_residue( false );
		// residue between built MS2 hairpin and the base pair step.
		add_mover.apply( pose, StepWiseMove( 86, Attachment(85,BOND_TO_PREVIOUS), ADD ) );
		// join the two poses
		add_mover.apply( pose, StepWiseMove( 87, Attachment(86,BOND_TO_PREVIOUS), ADD ) );

		// First residues 1-83 are protein. Rest are RNA:
		TS_ASSERT_EQUALS( pose.sequence(), "FANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPIFANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELTacaugaggaucacccagu" );

		utility::vector1< Size > working_moving_res_list = make_vector1( 86 );
		StepWiseWorkingParametersOP working_parameters = setup_working_parameters_for_swa( working_moving_res_list, pose, 0 /* no native pose */,  vector1< Size >() /*bridge_res*/, vector1< Size >()/*working_minimize_res*/  );

		vector1< PoseOP > pose_list;
		pose_list.push_back( pose.clone() );
		StepWiseModelerOptionsOP options( new StepWiseModelerOptions );
		TS_ASSERT_EQUALS( options->num_pose_minimize(), 0 ); // seek default

		StepWisePacker stepwise_packer( working_moving_res_list );
		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_res_level_energy.wts" ); // RNA/protein.
		stepwise_packer.set_scorefxn( scorefxn );
		TS_ASSERT( !stepwise_packer.pack_all_side_chains() );
		stepwise_packer.figure_out_neighbors( pose );
		TR << "working_pack_res: " << make_tag_with_dashes( stepwise_packer.working_pack_res() ) << std::endl;

		{
			// in enumeration mode, num_pose_minimize should be 0 --> minimize all
			TS_ASSERT( !options->choose_random() );
			StepWiseMinimizer stepwise_minimizer( pose_list, working_parameters, options, 0 /*dummy scorefxn*/  );
			stepwise_minimizer.set_working_pack_res( stepwise_packer.working_pack_res() );
			TS_ASSERT_EQUALS( stepwise_minimizer.num_pose_minimize(), 0 ); // default for RNA
			TS_ASSERT_EQUALS( stepwise_minimizer.working_minimize_res(), make_vector1( Size(86) ) ); // default for RNA

			TS_ASSERT( stepwise_minimizer.working_pack_res().has_value( Size(86) ) );  // 2'-OH for nearby residues
			TS_ASSERT( stepwise_minimizer.working_pack_res().has_value( Size(58) ) );  // protein side-chain for nearby residues
			TS_ASSERT_LESS_THAN( stepwise_minimizer.working_pack_res().size(), pose.total_residue() );  // but do not pack all residues.
		}

		{
			// in choose_random mode, num_pose_minimize should be 1 (for RNA)
			options->set_choose_random( true );
			TS_ASSERT( options->choose_random() );
			StepWiseMinimizer stepwise_minimizer( pose_list, working_parameters, options, 0 /*dummy scorefxn*/  );
			TS_ASSERT_EQUALS( stepwise_minimizer.num_pose_minimize(), 1 ); // default for RNA
			TS_ASSERT_EQUALS( stepwise_minimizer.working_minimize_res(), make_vector1( Size(86) ) ); // default for RNA
			TS_ASSERT( options->choose_random() );
		}

		{
			// in choose_random mode, num_pose_minimize should be 5 if moving at least one protein residue
			StepWiseWorkingParametersOP working_parameters_minimize_protein = setup_working_parameters_for_swa( vector1<Size>() /*no sampled res*/, pose, 0 /* no native pose */,  vector1< Size >() /*bridge_res*/, make_vector1(1,2) /*minimize protein res*/  );

			TS_ASSERT( options->choose_random() );
			StepWiseMinimizer stepwise_minimizer( pose_list, working_parameters_minimize_protein, options, 0 /*dummy scorefxn*/  );
			TS_ASSERT_EQUALS( stepwise_minimizer.num_pose_minimize(), 5 ); // default for RNA
			TS_ASSERT_EQUALS( stepwise_minimizer.working_minimize_res(), make_vector1( Size( 1 ), Size( 2 ), Size( 86 ) ) );
			TS_ASSERT( options->choose_random() );
		}

	}



};



