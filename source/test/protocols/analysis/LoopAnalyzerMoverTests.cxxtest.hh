// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/analysis/LoopAnalyzerMoverTests.cxxtest.hh
/// @brief  tests for LoopAnalyzerMover: that a loop returns values in a particular range, and that it behaves properly with loop definitions near termini
/// @author Steven Lewis (smlewi@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/analysis/LoopAnalyzerMover.hh>

#include <protocols/loops/Loops.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>

static THREAD_LOCAL basic::Tracer TR("LoopAnalyzerMoverTests");

class LoopAnalyzerMoverTests : public CxxTest::TestSuite {
	//Define Variables

	core::pose::PoseOP poseop;
	core::Real const TOLERATED_ERROR = 1e-4;

public:

	void setUp(){
		core_init();

		if ( !poseop ) {
			poseop = create_trpcage_ideal_poseop();
		}

	}

	void tearDown(){

	}

	//pick two groups of residues from the middle of the Pose, call loops, and check that LAM runs
	void test_simple_loop(){

		//set up loops
		protocols::loops::Loops middle_loops;
		//Pose goes from 1 to 20
		middle_loops.add_loop(3, 5, 3, 0, 0);
		middle_loops.add_loop(13, 15, 13, 0, 0);

		//construct and run LAM.  false is for tracer reporting (unwanted here)
		protocols::analysis::LoopAnalyzerMover LAM(middle_loops, false);
		LAM.apply(*poseop);

		//Would be awesome to TS_ASSERT that the pose has not changed, but there's no good way to do that.

		//Would be nice to test the nicely formatted output LAM dumps to PDBs, but that's tied up in JD2 and beyond the scope of this test

		//find out what the results are block
		if ( false ) {

			utility::vector1< core::Real > const chainbreak_scores_test(LAM.get_chainbreak_scores());
			TR << "cb vector size " << chainbreak_scores_test.size() << std::endl;

			for ( core::Size index(1); index <= chainbreak_scores_test.size(); ++index ) {
				TR << "cb vector entry " << index << " " <<  chainbreak_scores_test[index] << std::endl;
			}

			TR << " total_score "    << LAM.get_total_score() << std::endl;
			TR << " max_rama "       << LAM.get_max_rama() << std::endl;
			TR << " max_chainbreak " << LAM.get_max_omega() << std::endl;
			TR << " max_omega "      << LAM.get_max_pbond() << std::endl;
			TR << " max_pbond "      << LAM.get_max_chainbreak() << std::endl;
		}

		//hardcoded results - these are more "integration test" style than carefully manually calculated

		//size is 10.  3 each in two loops, plus one residue on either side of each loop
		utility::vector1< core::Real > const chainbreak_scores = {0.11566, 0.370071, 0.000137218, 0.0978375, 0.122621, 0.6909, 0.319962, 0.998043, 0.255284, 0.190627};

		core::Real const total_score( -30.2765 );
		core::Real const max_rama( 1.1919 );
		core::Real const max_chainbreak( 0.998043 );
		core::Real const max_omega( 1.8394 );
		core::Real const max_pbond( -3.42852 );

		//do comparisons
		utility::vector1< core::Real > const chainbreak_scores_test(LAM.get_chainbreak_scores());

		TS_ASSERT_EQUALS(chainbreak_scores.size(), chainbreak_scores_test.size());
		for ( core::Size index(1); index <= chainbreak_scores.size(); ++index ) {
			TS_ASSERT_DELTA(chainbreak_scores[index], chainbreak_scores_test[index], TOLERATED_ERROR);
		}

		TS_ASSERT_DELTA( total_score,    LAM.get_total_score(), TOLERATED_ERROR);
		TS_ASSERT_DELTA( max_rama,       LAM.get_max_rama(), TOLERATED_ERROR);
		TS_ASSERT_DELTA( max_chainbreak, LAM.get_max_chainbreak(), TOLERATED_ERROR);
		TS_ASSERT_DELTA( max_omega,      LAM.get_max_omega(), TOLERATED_ERROR);
		TS_ASSERT_DELTA( max_pbond,      LAM.get_max_pbond(), TOLERATED_ERROR);
	}

	///@details LAM automatically extends the analysis region by 1 on either side; this needs to NOT happen when adjacent to C termini.  Primary test here is just 'do not crash'
	void test_Cterm_adjacent(){

		//set up loops
		protocols::loops::Loops cterm_loops;
		//Pose goes from 1 to 20
		cterm_loops.add_loop(17, 19, 17, 0, 0);

		//construct and run LAM.  false is for tracer reporting (unwanted here)
		protocols::analysis::LoopAnalyzerMover LAM(cterm_loops, false);
		LAM.apply(*poseop);

		//just test that the chainbreak scores vector is loop_length+1 (not +2)
		TS_ASSERT_EQUALS(LAM.get_chainbreak_scores().size(), 4);

	}

	///@details LAM needs to refuse to test the C-terminal residue
	void test_internal_terminiCterm(){

		//set up loops
		protocols::loops::Loops cterm_loops;
		//Pose goes from 1 to 20
		cterm_loops.add_loop(18, 20, 18, 0, 0);

		//construct and run LAM.  false is for tracer reporting (unwanted here)
		protocols::analysis::LoopAnalyzerMover LAM(cterm_loops, false);
		LAM.apply(*poseop);

		//just test that the chainbreak scores vector is (17, 18, 19)
		TS_ASSERT_EQUALS(LAM.get_chainbreak_scores().size(), 3);

	}


	///@details Nterm adjacent is "normal", because the peptide bond points to i+1 - we are just testing that the number of tested residues is as expected
	void test_Nterm_adjacent(){

		//set up loops
		protocols::loops::Loops term_loops;
		//Pose goes from 1 to 20
		term_loops.add_loop(2, 4, 4, 0, 0);

		//construct and run LAM.  false is for tracer reporting (unwanted here)
		protocols::analysis::LoopAnalyzerMover LAM(term_loops, false);
		LAM.apply(*poseop);

		//just test that the chainbreak scores vector is loop_length+2 - asymmetric with C-termini because peptide bond points to i+1 not i-1
		TS_ASSERT_EQUALS(LAM.get_chainbreak_scores().size(), 5);

	}

	//just test that the chainbreak scores vector is loop_length+1 - cannot test residue 0!
	void test_at_Nterm(){

		//set up loops
		protocols::loops::Loops term_loops;
		//Pose goes from 1 to 20
		term_loops.add_loop(1, 3, 3, 0, 0);

		//construct and run LAM.  false is for tracer reporting (unwanted here)
		protocols::analysis::LoopAnalyzerMover LAM(term_loops, false);
		LAM.apply(*poseop);

		//just test that the chainbreak scores vector is loop_length+1 - cannot test residue 0!
		TS_ASSERT_EQUALS(LAM.get_chainbreak_scores().size(), 4);

	}

	//test near termini that are not also pose boundaries
	void test_internal_termini(){

		core::pose::Pose double_pose(*poseop);
		//core::pose::Pose copy(double_pose); //I don't think append_pose_to_pose is safe if poses are same
		core::pose::append_pose_to_pose(double_pose, double_pose);
		//set up loops
		protocols::loops::Loops ok_loops;
		//Pose goes from 1 to 20, 21 to 40
		ok_loops.add_loop(3, 5, 3, 0, 0);
		ok_loops.add_loop(13, 15, 13, 0, 0);
		ok_loops.add_loop(23, 25, 23, 0, 0);
		ok_loops.add_loop(33, 35, 33, 0, 0);

		//construct and run LAM.  false is for tracer reporting (unwanted here)
		protocols::analysis::LoopAnalyzerMover LAMok(ok_loops, false);
		LAMok.apply(double_pose);

		//just test that the chainbreak scores vector is num_loops*(loop_length+2) 4 * (3+2)
		TS_ASSERT_EQUALS(LAMok.get_chainbreak_scores().size(), 20);
		//////////////////////////////////////////////////////////////////////////////////
		//now test the 4 near-terminal loops
		//////////////////////////////////////////////////////////////////////////////////
		protocols::loops::Loops nearterm_loops1, nearterm_loops2,nearterm_loops3,nearterm_loops4;
		//Pose goes from 1 to 20, 21 to 40
		nearterm_loops1.add_loop(2, 4, 3, 0, 0); // area 5 residues
		nearterm_loops2.add_loop(17, 19, 17, 0, 0); // 4 residues
		nearterm_loops3.add_loop(22, 24, 23, 0, 0); // 5 residues
		nearterm_loops4.add_loop(37, 39, 37, 0, 0); // 4 residues

		protocols::analysis::LoopAnalyzerMover LAM_nearterm1(nearterm_loops1, false);
		protocols::analysis::LoopAnalyzerMover LAM_nearterm2(nearterm_loops2, false);
		protocols::analysis::LoopAnalyzerMover LAM_nearterm3(nearterm_loops3, false);
		protocols::analysis::LoopAnalyzerMover LAM_nearterm4(nearterm_loops4, false);

		//for these 4 near-terminal loops, check that they get the right number of residues analyzed
		LAM_nearterm1.apply(double_pose);
		TS_ASSERT_EQUALS(LAM_nearterm1.get_chainbreak_scores().size(), 5);
		LAM_nearterm2.apply(double_pose);
		TS_ASSERT_EQUALS(LAM_nearterm2.get_chainbreak_scores().size(), 4);
		LAM_nearterm3.apply(double_pose);
		TS_ASSERT_EQUALS(LAM_nearterm3.get_chainbreak_scores().size(), 5);
		LAM_nearterm4.apply(double_pose);
		TS_ASSERT_EQUALS(LAM_nearterm4.get_chainbreak_scores().size(), 4);

		//////////////////////////////////////////////////////////////////////////////////
		//now test the 4 at-terminal loops
		//////////////////////////////////////////////////////////////////////////////////
		protocols::loops::Loops atterm_loops1, atterm_loops2,atterm_loops3,atterm_loops4;
		//Pose goes from 1 to 20, 21 to 40
		atterm_loops1.add_loop(1, 3, 3, 0, 0); // area 4 residues
		atterm_loops2.add_loop(18, 20, 18, 0, 0); // 3 residues
		atterm_loops3.add_loop(21, 23, 23, 0, 0); // 4 residues
		atterm_loops4.add_loop(38, 40, 38, 0, 0); // 3 residues

		protocols::analysis::LoopAnalyzerMover LAM_atterm1(atterm_loops1, false);
		protocols::analysis::LoopAnalyzerMover LAM_atterm2(atterm_loops2, false);
		protocols::analysis::LoopAnalyzerMover LAM_atterm3(atterm_loops3, false);
		protocols::analysis::LoopAnalyzerMover LAM_atterm4(atterm_loops4, false);

		//for these 4 near-terminal loops, check that they get the right number of residues analyzed
		LAM_atterm1.apply(double_pose);
		TS_ASSERT_EQUALS(LAM_atterm1.get_chainbreak_scores().size(), 4);
		LAM_atterm2.apply(double_pose);
		TS_ASSERT_EQUALS(LAM_atterm2.get_chainbreak_scores().size(), 3);
		LAM_atterm3.apply(double_pose);
		TS_ASSERT_EQUALS(LAM_atterm3.get_chainbreak_scores().size(), 4);
		LAM_atterm4.apply(double_pose);
		TS_ASSERT_EQUALS(LAM_atterm4.get_chainbreak_scores().size(), 3);

		////////////////////////////////////////////////////////////////////////////
		//try throwing some deliberately invalid loops
		////////////////////////////////////////////////////////////////////////////
		protocols::loops::Loops bad_loops1, bad_loops2;
		//Pose goes from 1 to 20, 21 to 40
		bad_loops1.add_loop(17, 23, 23, 0, 0); // crosses chainbreak
		bad_loops2.add_loop(48, 200, 180, 0, 0); // out of bounds

		protocols::analysis::LoopAnalyzerMover LAM_bad1(bad_loops1, false);
		protocols::analysis::LoopAnalyzerMover LAM_bad2(bad_loops2, false);

		try {
			set_throw_on_next_assertion_failure();
			LAM_bad1.apply(double_pose);
			TS_ASSERT( false ); //should have thrown
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string const expected_error_message = "Loops::Error, loop crosses chain boundary";
			std::size_t const found = e.msg().find(expected_error_message); //look for the substring
			TS_ASSERT_DIFFERS( found, std::string::npos ); //assert we did not fail to find
		}
		try {
			set_throw_on_next_assertion_failure();
			LAM_bad2.apply(double_pose);
			TS_ASSERT( false ); //should have thrown
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string const expected_error_message = "LoopRebuild::ERROR Loop definition out of boundary \n";
			std::size_t const found = e.msg().find(expected_error_message); //look for the substring
			TS_ASSERT_DIFFERS( found, std::string::npos ); //assert we did not fail to find
		}
	}

};
