// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/stepwise/monte_carlo/mover//StepWiseMoveSelectorTest.cxxtest.hh
/// @brief  test StepWiseMoveSelector
/// @author Rhiju Das (rhiju@stanford.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project headers
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.hh>
#include <protocols/stepwise/monte_carlo/mover/options/StepWiseMoveSelectorOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>

// Protocol Headers
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.hh>

#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("StepWiseMoveSelectorTest");


class StepWiseMoveSelectorTest : public CxxTest::TestSuite {
	//Define Variables

public:

	core::pose::PoseCOP pose;
	protocols::stepwise::monte_carlo::submotif::SubMotifLibraryCOP submotif_library;
	protocols::stepwise::monte_carlo::mover::StepWiseMoveSelectorOP stepwise_move_selector;

	void setUp(){

		core_init_with_additional_options( "-fasta protocols/stepwise/monte_carlo/mover/t_loop.fasta -s protocols/stepwise/monte_carlo/mover/up_to_ua_handle_t_loop_3l0u_RNA.pdb" );

		using namespace core::chemical;
		using namespace core::pose;
		using namespace core::import_pose;
		using namespace protocols::stepwise::monte_carlo::submotif;
		using namespace protocols::stepwise::monte_carlo::mover;
		using namespace protocols::stepwise::monte_carlo::mover::options;
		using namespace utility;

		ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		pose = initialize_pose_and_other_poses_from_command_line( rsd_set );

		submotif_library = SubMotifLibraryOP( new SubMotifLibrary( rsd_set ) );

		StepWiseMoveSelectorOptionsOP options( new StepWiseMoveSelectorOptions );
		options->set_submotif_frequency( 0.2 );
		options->set_from_scratch_frequency( 0.1 );
		options->set_docking_frequency( 0.0 );

		stepwise_move_selector = StepWiseMoveSelectorOP( new StepWiseMoveSelector( options ) );
		stepwise_move_selector->set_submotif_library( submotif_library );
	}

	void tearDown(){
	}

	void test_submotif_split() {

		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace protocols::stepwise::monte_carlo::mover;
		using namespace protocols::stepwise::monte_carlo::mover::options;
		using namespace utility::tools;

		stepwise_move_selector->figure_out_all_possible_moves( *pose );
		TS_ASSERT( stepwise_move_selector->swa_moves().has_value( StepWiseMove( 5, Attachment( 4, BOND_TO_PREVIOUS ), ADD ) ) );

		StepWiseMove uturn_addition_move1( make_vector1( 5,6,7 ), make_vector1( Attachment( 4, BOND_TO_PREVIOUS ) ), ADD_SUBMOTIF, "uturns/uturn_t_loop_3l0u_RNA.pdb" );
		TS_ASSERT( stepwise_move_selector->swa_moves().has_value( uturn_addition_move1 ) );

		StepWiseMove uturn_addition_move2( make_vector1( 5,6,7 ), make_vector1( Attachment( 8, BOND_TO_NEXT ) ), ADD_SUBMOTIF, "uturns/uturn_t_loop_3l0u_RNA.pdb" );
		TS_ASSERT( stepwise_move_selector->swa_moves().has_value( uturn_addition_move2 ) );


		// add U-turn submotif
		PoseOP pose_add_submotif = pose->clone();
		StepWiseMove swa_move( uturn_addition_move1 );
		nonconst_full_model_info( *pose_add_submotif ).add_other_pose(
			submotif_library->create_new_submotif( swa_move.move_element(),
			swa_move.submotif_tag(), *pose_add_submotif ) );

		AddMover add_mover;
		add_mover.set_start_added_residue_in_aform( true );
		add_mover.set_presample_added_residue( false );
		add_mover.apply( *pose_add_submotif, swa_move );
		stepwise_move_selector->figure_out_all_possible_moves( *pose_add_submotif );

		StepWiseMove reverse_submotif_move = stepwise_move_selector->reverse_move( swa_move, *pose, *pose_add_submotif );  TS_ASSERT_EQUALS( reverse_submotif_move, StepWiseMove( make_vector1( 5,6,7 ), Attachment( 4, BOND_TO_PREVIOUS ), DELETE ) );

		StepWiseMove split_submotif_move( make_vector1( 6,7 ), Attachment( 5, BOND_TO_PREVIOUS ), DELETE );

		TS_ASSERT( stepwise_move_selector->swa_moves().has_value( reverse_submotif_move ) );
		TS_ASSERT( !stepwise_move_selector->swa_moves().has_value( split_submotif_move ) );

		StepWiseMoveSelectorOptionsOP options = stepwise_move_selector->options()->clone();
		options->set_allow_submotif_split( true );
		stepwise_move_selector->set_options( options );
		stepwise_move_selector->figure_out_all_possible_moves( *pose_add_submotif );
		TS_ASSERT( stepwise_move_selector->swa_moves().has_value( reverse_submotif_move ) );
		TS_ASSERT( stepwise_move_selector->swa_moves().has_value( split_submotif_move ) );


	}

	void test_reverse_of_FROM_SCRATCH_move(){

		using namespace core::pose;
		using namespace protocols::stepwise::monte_carlo::mover;
		using namespace utility::tools;

		stepwise_move_selector->figure_out_all_possible_moves( *pose );
		StepWiseMove swa_move( make_vector1(5,6), Attachments(), FROM_SCRATCH );
		TS_ASSERT( stepwise_move_selector->swa_moves().has_value( swa_move ) );

		// do FROM_SCRATCH
		PoseOP pose_do_from_scratch = pose->clone();
		FromScratchMover from_scratch_mover;
		from_scratch_mover.apply( *pose_do_from_scratch, swa_move.move_element() );

		// // check for reverse
		StepWiseMove reverse_from_scratch_move = stepwise_move_selector->reverse_move( swa_move, *pose, *pose_do_from_scratch );
		TR << reverse_from_scratch_move << std::endl;
		TS_ASSERT_EQUALS( reverse_from_scratch_move, StepWiseMove( 6, Attachment( 5, BOND_TO_PREVIOUS ), DELETE ) );
		stepwise_move_selector->figure_out_all_possible_moves( *pose_do_from_scratch );
		TS_ASSERT( stepwise_move_selector->swa_moves().has_value( reverse_from_scratch_move ) );


	}

};



