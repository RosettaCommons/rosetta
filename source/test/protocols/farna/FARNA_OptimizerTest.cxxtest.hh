// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/farna//FARNA_OptimizerTest.cxxtest.hh
/// @brief  test FARNA_Optimizer, wrapper for calling FARNA inside stepwise
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
#include <core/id/NamedAtomID.hh>

// Protocol Headers
#include <protocols/farna/FARNA_Optimizer.hh>
#include <protocols/farna/RNA_FragmentMonteCarlo.hh>
#include <protocols/farna/libraries/RNA_ChunkLibrary.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>

#include <utility/string_util.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("FARNA_OptimizerTest");


class FARNA_OptimizerTest : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-s protocols/farna/srl_fixed_START1_1q9a_RNA.pdb -terminal_res  A:651 A:669  -extra_min_res  A:652 A:658 A:663 A:668  -jump_res  A:651 A:669 A:652 A:668 A:658 A:663 A:659 A:662  -cutpoint_closed  A:651 A:658-659  -fasta protocols/farna/farna_rebuild.fasta" );
	}

	void tearDown(){

	}

	void test_sarcin_ricin_loop_buildup(){
		/////////////////////////////////////////////////////
		//
		//    _________________ Jump___
		//   |                         |
		//   1x (2)                (8)x9x 10
		//   |   |                  |  |
		//  19 (18)                (13)12 11
		//
		//  x mark cutpoints_closed
		//
		//  all residues come from starter PDB, but ones in parentheses can
		//  be moved during optimization (-extra_min_res)
		//
		// Numbering above is 'full_model' numbering.
		// In starter pose, numbers will be:
		//
		//    _________________ Jump___
		//   |                         |
		//   1x (2)                (3)x4x 5
		//   |   |                  |  |
		//  10  (9)                (8) 7 6
		//
		/////////////////////////////////////////////////////
		// Following is 'standard' setup from stepwise.cc

		using namespace core;
		using namespace core::chemical;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		using namespace core::id;
		using namespace protocols::farna;
		using namespace protocols::farna::libraries;
		using namespace protocols::toolbox;
		using namespace protocols::stepwise::setup;

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		PoseOP pose_op = initialize_pose_and_other_poses_from_command_line( rsd_set );
		pose::Pose & pose = *pose_op;
		TS_ASSERT( pose.fold_tree().is_cutpoint( 1 ) );
		TS_ASSERT( pose.fold_tree().is_cutpoint( 2 ) );
		TS_ASSERT( pose.fold_tree().is_cutpoint( 3 ) );
		TS_ASSERT( pose.fold_tree().is_cutpoint( 4 ) );
		TS_ASSERT( !pose.residue_type(1).has_variant_type( CUTPOINT_UPPER ) );
		TS_ASSERT(  pose.residue_type(1).has_variant_type( CUTPOINT_LOWER ) );
		TS_ASSERT(  pose.residue_type(2).has_variant_type( CUTPOINT_UPPER ) );
		TS_ASSERT( !pose.residue_type(2).has_variant_type( CUTPOINT_LOWER ) );
		TS_ASSERT( !pose.residue_type(3).has_variant_type( CUTPOINT_UPPER ) );
		TS_ASSERT(  pose.residue_type(3).has_variant_type( CUTPOINT_LOWER ) );
		TS_ASSERT(  pose.residue_type(4).has_variant_type( CUTPOINT_UPPER ) );
		TS_ASSERT(  pose.residue_type(4).has_variant_type( CUTPOINT_LOWER ) );
		TS_ASSERT(  pose.residue_type(5).has_variant_type( CUTPOINT_UPPER ) );
		TS_ASSERT( !pose.residue_type(5).has_variant_type( CUTPOINT_LOWER ) );
		TS_ASSERT( const_full_model_info( pose ).res_list() == utility::vector1< Size >( utility::get_resnum_and_chain( "1-2 8-13 18-19" ).first ) );

		Pose start_pose = pose;
		{

			utility::vector1< PoseOP > pose_list;
			pose_list.push_back( pose_op->clone() );

			ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_lores_for_stepwise.wts" );
			FARNA_Optimizer farna_optimizer( pose_list, lores_scorefxn, 5 /* cycles */ );
			farna_optimizer.apply( pose );

			// there should only be changes at extra_min_res.
			AtomLevelDomainMap const & atom_level_domain_map = *(farna_optimizer.rna_fragment_monte_carlo()->atom_level_domain_map());
			atom_level_domain_map.show( TR );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 1 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 2 ), pose ), ROSETTA_LIBRARY_DOMAIN );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 3 ), pose ), ROSETTA_LIBRARY_DOMAIN );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 4 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 5 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 6 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 7 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 8 ), pose ), ROSETTA_LIBRARY_DOMAIN );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 9 ), pose ), ROSETTA_LIBRARY_DOMAIN );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",10 ), pose ), 1 );

			TS_ASSERT_DELTA( pose.residue( 1 ).xyz( " C1'" ), start_pose.residue( 1 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 4 ).xyz( " C1'" ), start_pose.residue( 4 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 5 ).xyz( " C1'" ), start_pose.residue( 5 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 6 ).xyz( " C1'" ), start_pose.residue( 6 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 7 ).xyz( " C1'" ), start_pose.residue( 7 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue(10 ).xyz( " C1'" ), start_pose.residue(10 ).xyz( " C1'" ), 0.001 );

			TS_ASSERT_DELTA( pose.residue( 4 ).xyz( "OVL1" ), start_pose.residue( 4 ).xyz( "OVL1" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 4 ).xyz( "OVL2" ), start_pose.residue( 4 ).xyz( "OVL2" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 5 ).xyz( "OVU1" ), start_pose.residue( 5 ).xyz( "OVU1" ), 0.001 );
		}

		////////////////////////////////////////////////////////////////////////
		// Now let's fill in the rest of the RNA
		//
		//    _________________ Jump___
		//   |                         |
		//   1x (2)--3---4--5  6--7---(8)x9x 10
		//   |   |                         |  |
		//  19 (18) 17--16-----15-14-(13)12 11
		////////////////////////////////////////////////////////////////////////
		using namespace protocols::stepwise::monte_carlo::mover;
		AddMover add_mover;
		add_mover.set_start_added_residue_in_aform( true );
		add_mover.set_presample_added_residue( false );
		add_mover.apply( pose, StepWiseMove( 3, Attachment(2,BOND_TO_PREVIOUS), ADD ) );
		add_mover.apply( pose, StepWiseMove( 4, Attachment(3,BOND_TO_PREVIOUS), ADD ) );
		add_mover.apply( pose, StepWiseMove( 5, Attachment(4,BOND_TO_PREVIOUS), ADD ) );
		add_mover.apply( pose, StepWiseMove( 7, Attachment(8,BOND_TO_NEXT    ), ADD ) );
		add_mover.apply( pose, StepWiseMove( 6, Attachment(7,BOND_TO_NEXT    ), ADD ) );
		add_mover.apply( pose, StepWiseMove( 14, Attachment(13,BOND_TO_PREVIOUS ), ADD ) );
		add_mover.apply( pose, StepWiseMove( 15, Attachment(14,BOND_TO_PREVIOUS ), ADD ) );
		add_mover.apply( pose, StepWiseMove( 16, Attachment(15,BOND_TO_PREVIOUS ), ADD ) );
		add_mover.apply( pose, StepWiseMove( 17, Attachment(16,BOND_TO_PREVIOUS ), ADD ) );


		{
			utility::vector1< PoseOP > pose_list;
			pose_list.push_back( pose.clone() );

			ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_lores_for_stepwise.wts" );
			FARNA_Optimizer farna_optimizer( pose_list, lores_scorefxn, 5 /* cycles */ );
			farna_optimizer.apply( pose );

			// there should only be changes at extra_min_res.
			AtomLevelDomainMap const & atom_level_domain_map = *(farna_optimizer.rna_fragment_monte_carlo()->atom_level_domain_map());
			atom_level_domain_map.show( TR );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 1 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 2 ), pose ), ROSETTA_LIBRARY_DOMAIN );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 3 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 4 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 5 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 6 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 7 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 8 ), pose ), ROSETTA_LIBRARY_DOMAIN );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'", 9 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",10 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",11 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",12 ), pose ), 1 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",13 ), pose ), ROSETTA_LIBRARY_DOMAIN );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",14 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",15 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",16 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",17 ), pose ), 0 );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",18 ), pose ), ROSETTA_LIBRARY_DOMAIN );
			TS_ASSERT_EQUALS( atom_level_domain_map.get_domain( NamedAtomID( " C1'",19 ), pose ), 1 );

			TS_ASSERT_DELTA( pose.residue( 1 ).xyz( " C1'" ), start_pose.residue( 1 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 9 ).xyz( " C1'" ), start_pose.residue( 4 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue(10 ).xyz( " C1'" ), start_pose.residue( 5 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue(11 ).xyz( " C1'" ), start_pose.residue( 6 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue(12 ).xyz( " C1'" ), start_pose.residue( 7 ).xyz( " C1'" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue(19 ).xyz( " C1'" ), start_pose.residue(10 ).xyz( " C1'" ), 0.001 );

			TS_ASSERT_DELTA( pose.residue( 9 ).xyz( "OVL1" ), start_pose.residue( 4 ).xyz( "OVL1" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 9 ).xyz( "OVL2" ), start_pose.residue( 4 ).xyz( "OVL2" ), 0.001 );
			TS_ASSERT_DELTA( pose.residue( 10 ).xyz( "OVU1" ), start_pose.residue( 5 ).xyz( "OVU1" ), 0.001 );
		}

	}



};



