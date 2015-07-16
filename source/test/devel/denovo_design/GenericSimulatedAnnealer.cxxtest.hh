// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/denovo_design/GenericSimulatedAnnealer.cxxtest.hh
/// @brief  test suite for devel::denovo_design::GenericSimulatedAnnealer
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <devel/denovo_design/GenericSimulatedAnnealer.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/pointer/owning_ptr.hh>

/// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/simple_filters/ResidueCountFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR("devel.denovo_design.GenericSimulatedAnnealer.cxxtest");

// --------------- Test Class --------------- //
class GenericSimulatedAnnealerTests : public CxxTest::TestSuite {
	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn;

	// filter that simply returns the number of TRP residues
	utility::pointer::shared_ptr< protocols::simple_filters::ResidueCountFilter > num_trp;
	utility::pointer::shared_ptr< protocols::simple_filters::ScoreTypeFilter > score_filter;

	// mutate residues to trp and ala
	protocols::simple_moves::MutateResidueOP mut_to_trp, mut_to_trp2;
	protocols::simple_moves::MutateResidueOP mut_to_ala, mut_to_ala2;
	protocols::simple_moves::MutateResidueOP mut_to_phe;

	// original pose read from disk
	core::pose::Pose input_pose;

public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSet & residue_set = ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD );
		if ( !residue_set.has_name("D2I") )
			params_files.push_back("devel/denovo_design/D2I.params");
		residue_set.read_files(params_files);

		// initialize common filters/movers/scorefxns
		scorefxn = core::scoring::get_score_function( true );

		num_trp = utility::pointer::shared_ptr< protocols::simple_filters::ResidueCountFilter >(
				new protocols::simple_filters::ResidueCountFilter() );
		utility::vector1< std::string > res_to_count;
		res_to_count.push_back( "TRP" );
		num_trp->res_types( res_to_count );

		// for test purposes, we will mutate residues 25 and 52 to TRP
		mut_to_trp = protocols::simple_moves::MutateResidueOP( new protocols::simple_moves::MutateResidue( 25, "TRP" ) );
		mut_to_trp2 = protocols::simple_moves::MutateResidueOP( new protocols::simple_moves::MutateResidue( 52, "TRP" ) );
		// mutate residues 25 and 52 to ALA
		mut_to_ala = protocols::simple_moves::MutateResidueOP( new protocols::simple_moves::MutateResidue( 25, "ALA" ) );
		mut_to_ala2 = protocols::simple_moves::MutateResidueOP( new protocols::simple_moves::MutateResidue( 52, "ALA" ) );
		// mutate residue 25 to PHE
		mut_to_phe = protocols::simple_moves::MutateResidueOP( new protocols::simple_moves::MutateResidue( 25, "PHE" ) );

		// get the total score
		score_filter = utility::pointer::shared_ptr< protocols::simple_filters::ScoreTypeFilter >(
				new protocols::simple_filters::ScoreTypeFilter( scorefxn,
					core::scoring::score_type_from_name( "total_score" ),
					999999.9 ) );

		std::string const pdb_file( "devel/denovo_design/test_input.pdb" );
		core::import_pose::pose_from_pdb( input_pose, pdb_file );
		mut_to_trp->apply( input_pose );
		mut_to_trp2->apply( input_pose );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// test overall acceptance criteria code
	void test_acceptance() {
		TR << "Starting test_acceptance()" << std::endl;
		core::pose::Pose pose( input_pose );

		devel::denovo_design::GenericSimulatedAnnealer annealer;
		annealer.set_mover( mut_to_ala );
		annealer.add_filter( num_trp, true, 0.0, "low", true );
		annealer.set_maxtrials(2);
		annealer.apply( pose );
		TS_ASSERT( annealer.num_accepted_scores() == 3 );
		utility::vector1< utility::vector1< core::Real > > accepted_scores;
		for ( core::Size i=1; i<=annealer.num_accepted_scores(); ++i ) {
			accepted_scores.push_back( annealer.accepted_scores( i ) );
		}
		// Should be 3 (start pose), 2 (for 1st iteration), 2 (for 2nd iteration)
		TS_ASSERT_DELTA( accepted_scores[1][1], 3.000000000, 1e-9 );
		TS_ASSERT_DELTA( accepted_scores[2][1], 2.000000000, 1e-9 );
		TS_ASSERT_DELTA( accepted_scores[3][1], 2.000000000, 1e-9 );
		// Best score should be returned
		TS_ASSERT_DELTA( num_trp->report_sm( pose ), annealer.lowest_scores()[1], 1e-6 );

		// now test a rejection -- both runs of this mover should be rejected
		core::pose::Pose rej_pose( pose );
		annealer.set_mover( mut_to_trp );
		annealer.apply( rej_pose );
		TS_ASSERT( annealer.num_accepted_scores() == 1 );
		TS_ASSERT_DELTA( num_trp->report_sm( pose ), num_trp->report_sm( rej_pose ), 1e-6 );
		core::pose::Pose control_pose( pose );
		mut_to_trp->apply( control_pose );
		TS_ASSERT_DELTA( num_trp->report_sm( pose ) + 1, num_trp->report_sm( control_pose ), 1e-6 );
	}

	// check to see whether a file exists
	bool exists( std::string const & filename ) {
		std::ifstream f( filename.c_str() );
		bool const retval( f.good() );
		f.close();
		return retval;
	}

	// test the checkpointing code
	void test_checkpointing() {
		TR << "STARTING checkpointing test" << std::endl;

		// set up two genericsimulatedannealers
		std::string const checkpoint_file( "devel/denovo_design/mc" );
		devel::denovo_design::GenericSimulatedAnnealer annealer;
		annealer.checkpoint_file( checkpoint_file );
		annealer.set_mover( mut_to_ala );
		annealer.add_filter( num_trp, true, 0.0, "low", true );
		annealer.keep_checkpoint_file( true );
		TR << "Annealer filters: " << annealer.filters().size() << std::endl;
		annealer.set_maxtrials(1);

		devel::denovo_design::GenericSimulatedAnnealer annealer2;
		annealer2.checkpoint_file( checkpoint_file );
		annealer2.set_mover( mut_to_ala2 );
		annealer2.add_filter( num_trp, true, 0.0, "low", true );
		annealer2.keep_checkpoint_file( true );
		annealer2.set_maxtrials(1);

		core::pose::Pose pose( input_pose );

		// there should now be 3 TRP in the pose
		core::Real ntrp( num_trp->report_sm( pose ) );
		TR << "STARTING2" << std::endl;
		TR << "ntrp=" << ntrp << std::endl;
		TS_ASSERT_DELTA( ntrp, 3.000000000, 1e-9 );

		TR << "Clearing checkpoing data and starting annealer" << std::endl;
		annealer.remove_checkpoint_file();
		// run the first annealer
		annealer.apply( pose );
		TR << "Annealer finished" << std::endl;
		// there should be one accepted pose and the starting values
		TS_ASSERT( annealer.num_accepted_scores() == 2 );
		// the first trp should be mutated to ala
		ntrp = num_trp->report_sm( pose );
		TR << "ntrp after first annealer step =" << ntrp << std::endl;
		TS_ASSERT_DELTA( ntrp, 2.000000000, 1e-9 );
		// checkpoint files should be there.
		TS_ASSERT( exists( checkpoint_file ) );
		TS_ASSERT( exists(checkpoint_file + "_best.pdb") );
		TS_ASSERT( exists(checkpoint_file + "_last.pdb") );

		// best.pdb should contain 2 TRP, as should last.pdb
		core::pose::Pose test_pose;
		core::import_pose::pose_from_pdb( test_pose, checkpoint_file + "_best.pdb" );
		TS_ASSERT( num_trp->report_sm( test_pose ) == ntrp );
		core::import_pose::pose_from_pdb( test_pose, checkpoint_file + "_last.pdb" );
		TS_ASSERT( num_trp->report_sm( test_pose ) == ntrp );

		// when we run the second annealer, it should find the checkpoint file and not do anything
		annealer2.apply( pose );
		// there should be one accepted pose and the starting values in the accepted scores
		TS_ASSERT( annealer2.num_accepted_scores() == 2 );
		core::Real ntrp2 = num_trp->report_sm( pose );
		TR << "3ntrp=" << ntrp2 << std::endl;
		TS_ASSERT_DELTA( ntrp, ntrp2, 1e-9 );
		// and checkpoint files should still be there
		TS_ASSERT( exists(checkpoint_file) );
		TS_ASSERT( exists(checkpoint_file + "_best.pdb") );
		TS_ASSERT( exists(checkpoint_file + "_last.pdb") );

		// best.pdb should contain 2 TRP, as should last.pd
		core::import_pose::pose_from_pdb( test_pose, checkpoint_file + "_best.pdb" );
		TS_ASSERT( num_trp->report_sm( test_pose ) == ntrp2 );
		core::import_pose::pose_from_pdb( test_pose, checkpoint_file + "_last.pdb" );
		TS_ASSERT( num_trp->report_sm( test_pose ) == ntrp2 );

		// now we tell annealer2 to do two trials -- it should run once
		annealer2.set_maxtrials( 2 );
		annealer2.apply( pose );
		TS_ASSERT( annealer2.num_accepted_scores() == 3 );
		TS_ASSERT_DELTA( annealer2.accepted_scores(1)[1], 3.000000000, 1e-9 );
		TS_ASSERT_DELTA( annealer2.accepted_scores(2)[1], 2.000000000, 1e-9 );
		TS_ASSERT_DELTA( annealer2.accepted_scores(3)[1], 1.000000000, 1e-9 );
		core::Real ntrp3 = num_trp->report_sm( pose );
		TR << "after 2 trials ntrp=" << ntrp3 << std::endl;
		TS_ASSERT_DELTA( ntrp3, 1.000000000, 1e-9 );
		// checkpoint files should still be there
		TS_ASSERT( exists(checkpoint_file) );
		core::import_pose::pose_from_pdb( test_pose, checkpoint_file + "_best.pdb" );
		TS_ASSERT( num_trp->report_sm( test_pose ) == ntrp3 );
		core::import_pose::pose_from_pdb( test_pose, checkpoint_file + "_last.pdb" );
		TS_ASSERT( num_trp->report_sm( test_pose ) == ntrp3 );
		annealer.remove_checkpoint_file();
	}

	// test boltzmann acceptance function
	void test_boltzmann() {
		// copy input pose
		core::pose::Pose pose( input_pose );
		core::pose::Pose start_pose( pose );
		// now test with multiple filters
		devel::denovo_design::GenericSimulatedAnnealer annealer3;
		annealer3.set_maxtrials( 2 );
		annealer3.set_mover( mut_to_ala );
		// this one counts trp residues and scores
		annealer3.add_filter( num_trp, true, 0.01, "low", true );
		annealer3.add_filter( score_filter, true, 0.1, "low", false );
		annealer3.set_boltz_rank( true );
		// reset scoring and counters
		annealer3.reset( pose );
		// ensure the scoring counts are initialized properly
		TS_ASSERT( annealer3.num_accepted_scores() == 1 );
		TS_ASSERT( annealer3.accepted_scores( 1 ).size() == annealer3.filters().size() );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[2], score_filter->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[2], score_filter->report_sm( pose ), 1e-6 );

		// this mutation should help both score and num_trp
		// create a set of "random" numbers that aren't actually random to allow deterministic testing of boltzmann()
		utility::vector1< core::Real > randoms( annealer3.filters().size(), 1.0 );
		mut_to_ala->apply( pose );
		// should be accepted
		devel::denovo_design::TrialResult result( annealer3.boltzmann_result( pose, randoms ) );
		TS_ASSERT( result == devel::denovo_design::ACCEPTED );
		// ensure the scoring counts are tabulated properly
		TR << "number of accepted scores is now: " << annealer3.num_accepted_scores() << std::endl;
		TS_ASSERT( annealer3.num_accepted_scores() == 2 );
		TS_ASSERT( annealer3.accepted_scores( 2 ).size() == annealer3.filters().size() );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[2], score_filter->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[2], score_filter->report_sm( pose ), 1e-6 );

		// this one should be rejected -- it will have a worse score and all randoms are zero
		mut_to_phe->apply( pose );
		result = annealer3.boltzmann_result( pose, randoms );
		TS_ASSERT( result == devel::denovo_design::REJECTED );
		TS_ASSERT( annealer3.num_accepted_scores() == 2 );
		TS_ASSERT( annealer3.accepted_scores( 2 ).size() == annealer3.filters().size() );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[2], score_filter->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[2], score_filter->report_sm( pose ), 1e-6 );
		// residue 25 should be restored to ALA, not PHE
		TS_ASSERT( pose.residue( 25 ).name3() == "ALA" );

		// now set the random for the scorefilter to 0, which should accept everything for that filter
		randoms[2] = 0.0;
		// make the mutation back to TRP, which should be rejected based on ntrp filter
		mut_to_trp->apply( pose );
		result = annealer3.boltzmann_result( pose, randoms );
		TS_ASSERT( result == devel::denovo_design::REJECTED );
		TS_ASSERT( annealer3.num_accepted_scores() == 2 );
		TS_ASSERT( annealer3.accepted_scores( 2 ).size() == annealer3.filters().size() );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[2], score_filter->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[2], score_filter->report_sm( pose ), 1e-6 );
		TS_ASSERT( pose.residue( 25 ).name3() == "ALA" );

		// the current pose should be the best we will ever see
		core::pose::Pose best_pose( pose );
		// now set the other random to 0, which should accept everything
		randoms[1] = 0.0;
		// make the mutation back to TRP, which should be accepted based on temperature
		mut_to_trp->apply( pose );
		result = annealer3.boltzmann_result( pose, randoms );
		TS_ASSERT( result == devel::denovo_design::ACCEPTED );
		TS_ASSERT( annealer3.num_accepted_scores() == 3 );
		TS_ASSERT( annealer3.accepted_scores( 3 ).size() == annealer3.filters().size() );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[2], score_filter->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[1], num_trp->report_sm( best_pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[2], score_filter->report_sm( best_pose ), 1e-6 );
		TS_ASSERT( pose.residue( 25 ).name3() == "TRP" );
		TS_ASSERT( annealer3.lowest_score_pose()->residue( 25 ).name3() == "ALA" );

		// save the lowest, last_accepted, and the entire accepted_scores table
		utility::vector1< core::Real > saved_lowest_scores( annealer3.lowest_scores() );
		utility::vector1< core::Real > saved_last_accepted_scores( annealer3.last_accepted_scores() );
		utility::vector1< utility::vector1< core::Real > > saved_accepted_scores;
		for ( core::Size i=1; i<=annealer3.num_accepted_scores(); ++i ) {
			saved_accepted_scores.push_back( annealer3.accepted_scores( i ) );
		}

		// save the "best score" and "last accepted score"
		core::Real saved_lowest_score( annealer3.lowest_score() );
		core::Real saved_last_accepted_score( annealer3.last_accepted_score() );
		TR << "best=" << saved_lowest_score << "; last=" << saved_last_accepted_score << std::endl;

		// test temperature scaling function
		core::Real const temp_factor( 0.5 );
		annealer3.scale_temperatures( temp_factor );
		TR << "best=" << annealer3.lowest_score() << "; last=" << annealer3.last_accepted_score() << std::endl;
		// best, last accepted scores should be scaled
		TS_ASSERT_DELTA( annealer3.lowest_score() * temp_factor, saved_lowest_score, 1e-4 );
		TS_ASSERT_DELTA( annealer3.last_accepted_score() * temp_factor, saved_last_accepted_score, 1e-4 );
		// lowest scores and last accepted scores should be the same
		TS_ASSERT( saved_lowest_scores == annealer3.lowest_scores() );
		TS_ASSERT( saved_last_accepted_scores == annealer3.last_accepted_scores() );

		// now switch res 25 to phe which will result in a new acceptance but with different temps
		mut_to_phe->apply( pose );
		result = annealer3.boltzmann_result( pose, randoms );
		TS_ASSERT( result == devel::denovo_design::ACCEPTED );
		TS_ASSERT( annealer3.num_accepted_scores() == 4 );
		TS_ASSERT( annealer3.accepted_scores( 4 ).size() == annealer3.filters().size() );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[1], num_trp->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.last_accepted_scores()[2], score_filter->report_sm( pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[1], num_trp->report_sm( best_pose ), 1e-6 );
		TS_ASSERT_DELTA( annealer3.lowest_scores()[2], score_filter->report_sm( best_pose ), 1e-6 );
		TS_ASSERT( pose.residue( 25 ).name3() == "PHE" );
		TS_ASSERT( annealer3.lowest_score_pose()->residue( 25 ).name3() == "ALA" );

		// now when we return the best pose, it should be the best one with the ALA
		annealer3.recover_low( pose );
		TS_ASSERT( pose.residue( 25 ).name3() == "ALA" );
	}

};
