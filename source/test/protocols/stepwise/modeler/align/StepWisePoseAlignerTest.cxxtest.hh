// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/stepwise/modeler/align//StepWisePoseAlignerTest.cxxtest.hh
/// @brief  test StepWisePoseAligner
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
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.hh>
#include <protocols/stepwise/monte_carlo/util.hh>

#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("StepWisePoseAlignerTest");


using namespace core::pose;

class StepWisePoseAlignerTest : public CxxTest::TestSuite {
	//Define Variables

public:

	PoseOP pose, native_pose, align_pose;

	void setUp(){
		core_init_with_additional_options( "-s protocols/stepwise/modeler/align/1_scaff_subset_stepwise_input_seq_1.pdb protocols/stepwise/modeler/align/2_scaff_subset_stepwise_input_seq_1.pdb  -fasta protocols/stepwise/modeler/align/seq_with_AAmismatch.fasta -align_pdb protocols/stepwise/modeler/align/scaff_subset_stepwise_input_seq_1.pdb -pack_missing_sidechains false " );
	}

	void tearDown(){
	}

	void test_superimpose_res_in_root_partition()
	{
		using namespace core::chemical;
		using namespace core::pose::full_model_info;
		using namespace protocols::stepwise::setup;
		using namespace protocols::stepwise::monte_carlo;
		using namespace protocols::stepwise::monte_carlo::mover;
		using namespace protocols::stepwise::modeler;
		using namespace protocols::stepwise::modeler::align;

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		PoseOP pose = initialize_pose_and_other_poses_from_command_line( rsd_set );

		PoseOP native_pose, align_pose;
		initialize_native_and_align_pose( native_pose, align_pose, rsd_set, pose );

		TR << get_all_res_list( *pose ) << std::endl;

		// a little RNA base pair step to build onto the end of the MS2 RNA hairpin.
		TS_ASSERT_EQUALS( pose->sequence(), "acgu" );

		FullModelInfo const & full_model_info = const_full_model_info( *pose );

		// the other pose is the protein with the tip of the RNA hairpin.
		TS_ASSERT_EQUALS( full_model_info.other_pose_list().size(), 1 );
		Pose const & other_pose = *full_model_info.other_pose_list()[ 1 ];
		TS_ASSERT_EQUALS( other_pose.sequence(),  "FANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPIFANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELTugaggaucaccca" );

		TS_ASSERT_EQUALS( native_pose, align_pose );
		TS_ASSERT_DIFFERS( align_pose, PoseOP( 0 ) );
		TS_ASSERT_EQUALS( align_pose->sequence(), "FANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPIFANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELTacugaggaucacccagu" );

		TS_ASSERT_EQUALS( full_model_info.full_sequence(), "FANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPIFANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELTacaugaggaucacccaagu" );


		// let's imagine joining these together. Then the pose
		// will have pieces inherited from two different input structures -- two possible
		// domains on which to align.
		AddMover add_mover;
		add_mover.set_start_added_residue_in_aform( true );
		add_mover.set_presample_added_residue( false );
		// residue between built MS2 hairpin and the base pair step.
		add_mover.apply( *pose, StepWiseMove( 86, Attachment(85,BOND_TO_PREVIOUS), ADD ) );
		// join the two poses
		add_mover.apply( *pose, StepWiseMove( 87, Attachment(86,BOND_TO_PREVIOUS), ADD ) );

		TS_ASSERT_EQUALS( pose->sequence(), "FANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPIFANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELTacaugaggaucacccagu" );


		TS_ASSERT_EQUALS( pose->fold_tree().root(), 84 );
		StepWisePoseAligner stepwise_pose_aligner( *align_pose );
		stepwise_pose_aligner.apply( *pose );
		TS_ASSERT_EQUALS( make_tag_with_dashes(stepwise_pose_aligner.superimpose_res_in_pose() ), "84-85 100-101" );

		// Let's move the root -- that's what would happen during stepwise modeling, which
		// prefers roots within bigger domains.
		FoldTree f = pose->fold_tree();
		f.reorder( 1 );
		pose->fold_tree( f );

		TS_ASSERT_EQUALS( pose->fold_tree().root(), 1 );
		stepwise_pose_aligner.apply( *pose );
		TS_ASSERT_EQUALS( make_tag_with_dashes(stepwise_pose_aligner.superimpose_res_in_pose() ), "1-83 87-99" );

	}

	void test_rna_dna_match()
	{
		using namespace core::chemical;
		using namespace core::import_pose;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace protocols::stepwise::setup;
		using namespace protocols::stepwise::modeler::align;
		using namespace utility::tools;

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		PoseOP pose = initialize_pose_and_other_poses_from_command_line( rsd_set );

		// the other pose is the protein with the tip of the RNA hairpin.
		FullModelInfo const & full_model_info = const_full_model_info( *pose );
		Pose & other_pose = *full_model_info.other_pose_list()[ 1 ];

		PoseOP dna_pose = get_pdb_with_full_model_info( "protocols/stepwise/modeler/align/scaff_subset_stepwise_input_seq_1_mutated_to_DNA.pdb", rsd_set );
		TS_ASSERT_EQUALS( other_pose.sequence(), "FANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPIFANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELTugaggaucaccca" );

		TS_ASSERT_EQUALS( dna_pose->sequence(),  "FANGVAEWKVTCSVRQSSAQNRKYTIKVEYLNMELTIPIFANGVAEWRSKVTCSVRQSSANRKYTIKVEVPAAWRSYLNMELTactgaggatcacccagt" );

		StepWisePoseAligner stepwise_pose_aligner( *dna_pose );
		stepwise_pose_aligner.set_user_defined_calc_rms_res( make_vector1(87,88,89,90) ); // RNA residues
		TS_ASSERT_EQUALS( other_pose.residue_type( 84 ).aa(), na_ura );
		TS_ASSERT_EQUALS( dna_pose->residue_type( 86 ).aa(),  na_thy );
		stepwise_pose_aligner.apply( other_pose );

		other_pose.dump_pdb( "OTHER_POSE.pdb" );
		dna_pose->dump_pdb( "DNA_POSE.pdb" );

		TS_ASSERT_LESS_THAN( stepwise_pose_aligner.rmsd(), 1.0e-4);

	}


};
