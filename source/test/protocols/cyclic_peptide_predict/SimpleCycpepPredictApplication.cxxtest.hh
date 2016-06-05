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
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI_JobResultsSummary.hh>
#include <protocols/cyclic_peptide_predict/util.hh>

// Other Rosetta libraries:
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/MutateResidue.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.cyclic_peptide_predict.simple_cycpep_predict.cxxtest.hh");


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

	/// @brief This is the simplest test in the world: run the protocol and see if it crashes.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_simple_cycpep_predict() {
		TR << "Running SimpleCycpepPredictApplicationTests::test_simple_cycpep_predict()." << std::endl;
		TR << "This simply runs the SimpleCycpepPredictApplication and checks for crashes." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;
		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication the_app;
		the_app.run();
	}

	/// @brief Test the sort algorithm used by the MPI version to sort job results by energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_quicksort() {
		TR << "Running SimpleCycpepPredictApplicationTests::test_quicksort()." << std::endl;
		TR << "This tests the quicksort algorithm used by the MPI version of the app." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;

		using namespace protocols::cyclic_peptide_predict;

		utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist;
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  1, -12.31, 1.23, 1) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  2, -14.87, 2.72, 3) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  3,  -5.12, 3.11, 5) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  4,   2.16, 4.52, 2) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  5, -13.99, 1.83, 0) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  6, -18.54, 0.07, 5) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  7, -16.12, 1.32, 2) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  8, -13.43, 0.32, 4) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  9, -10.63, 0.86, 3) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3, 10,  -8.12, 1.96, 2) ) );

		TR << "Before sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist.size(); i<=imax; ++i ) {
			summarylist[i]->show(TR); TR.flush();
		}
		TR << std::endl;

		sort_jobsummaries_list( summarylist, SORT_BY_ENERGIES );

		TR << "After sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist.size(); i<=imax; ++i ) {
			summarylist[i]->show(TR);  TR.flush();
			if ( i>1 ) TS_ASSERT( summarylist[i]->pose_energy() >= summarylist[i-1]->pose_energy() );
		}
		TR << std::endl;

	}

	/// @brief Test the sort algorithm used by the MPI version to sort job results by rmsd.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_quicksort_rmsd() {
		TR << "Running SimpleCycpepPredictApplicationTests::test_quicksort()." << std::endl;
		TR << "This tests the quicksort algorithm used by the MPI version of the app, sorting by rmsd." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;

		using namespace protocols::cyclic_peptide_predict;

		utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist;
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  1, -12.31, 1.23, 1) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  2, -14.87, 2.72, 4) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  3,  -5.12, 3.11, 3) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  4,   2.16, 4.52, 1) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  5, -13.99, 1.83, 7) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  6, -18.54, 0.07, 5) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  7, -16.12, 1.32, 6) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  8, -13.43, 0.32, 6) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  9, -10.63, 0.86, 4) ) );
		summarylist.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3, 10,  -8.12, 1.96, 3) ) );

		TR << "Before sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist.size(); i<=imax; ++i ) {
			summarylist[i]->show(TR); TR.flush();
		}
		TR << std::endl;

		sort_jobsummaries_list( summarylist, SORT_BY_RMSD );

		TR << "After sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist.size(); i<=imax; ++i ) {
			summarylist[i]->show(TR);  TR.flush();
			if ( i>1 ) TS_ASSERT( summarylist[i]->rmsd() >= summarylist[i-1]->rmsd() );
		}
		TR << std::endl;

	}

	/// @brief Test the sort algorithm used by the MPI version to sort job results by energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_mergesort() {
		TR << "Running SimpleCycpepPredictApplicationTests::test_mergesort()." << std::endl;
		TR << "This tests the merge sort algorithm used by the MPI version of the app." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;

		using namespace protocols::cyclic_peptide_predict;

		utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist1;
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  1, -12.31, 1.23, 7) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  2, -11.41, 2.72, 5) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  3,  -7.93, 3.11, 3) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  4,  -6.33, 4.52, 1) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  5,  -5.99, 1.83, 2) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  6,  -5.54, 0.07, 3) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  7,  -2.12, 1.32, 8) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  8,  -1.43, 0.32, 3) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  9,   1.63, 0.86, 1) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3, 10,   2.12, 1.96, 3) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3, 11,   2.13, 1.53, 2) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3, 12,   2.15, 1.54, 1) ) );

		utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist2;
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  1, -14.74, 7.23, 6) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  2, -13.86, 3.78, 7) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  3, -13.32, 4.51, 3) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  4, -10.28, 0.52, 5) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  5,  -9.97, 1.83, 2) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  6,  -8.11, 6.17, 6) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  7,  -6.10, 8.35, 2) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  8,  -3.09, 2.35, 1) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  9,  -0.06, 9.86, 2) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5, 10,   2.88, 4.46, 1) ) );

		TR << "List 1 before sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
			summarylist1[i]->show(TR); TR.flush();
		}
		TR << std::endl;

		TR << "List 2 before sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist2.size(); i<=imax; ++i ) {
			summarylist2[i]->show(TR); TR.flush();
		}
		TR << std::endl;

		mergesort_jobsummaries_list( summarylist1, summarylist2, SORT_BY_ENERGIES );

		TR << "After sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
			summarylist1[i]->show(TR);  TR.flush();
			if ( i>1 ) TS_ASSERT( summarylist1[i]->pose_energy() >= summarylist1[i-1]->pose_energy() );
		}
		TR << std::endl;

		TS_ASSERT_EQUALS( summarylist1.size(), 22);

	}

	/// @brief Test the sort algorithm used by the MPI version to sort job results by energy.
	/// @details This version tests the edge cases of empty lists.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_mergesort_edgecases() {
		TR << "Running SimpleCycpepPredictApplicationTests::test_mergesort()." << std::endl;
		TR << "This tests the merge sort algorithm used by the MPI version of the app.  This test tests the edge cases of empty lists." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;

		using namespace protocols::cyclic_peptide_predict;

		{ //2nd list empty
			TR << "Trying with second list empty..." << std::endl;
			utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist1;
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  1, -12.31, 1.23, 6) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  2, -11.41, 2.72, 7) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  3,  -7.93, 3.11, 4) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  4,  -6.33, 4.52, 7) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  5,  -5.99, 1.83, 5) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  6,  -5.54, 0.07, 3) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  7,  -2.12, 1.32, 5) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  8,  -1.43, 0.32, 2) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  9,   1.63, 0.86, 1) ) );
			summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3, 10,   2.12, 1.96, 3) ) );

			utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist2;

			TR << "List 1 before sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
				summarylist1[i]->show(TR); TR.flush();
			}
			TR << std::endl;

			TR << "List 2 before sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist2.size(); i<=imax; ++i ) {
				summarylist2[i]->show(TR); TR.flush();
			}
			TR << std::endl;

			mergesort_jobsummaries_list( summarylist1, summarylist2, SORT_BY_ENERGIES );

			TR << "After sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
				summarylist1[i]->show(TR);  TR.flush();
				if ( i>1 ) TS_ASSERT( summarylist1[i]->pose_energy() >= summarylist1[i-1]->pose_energy() );
			}
			TR << std::endl;
			TS_ASSERT_EQUALS(summarylist1.size(), 10);
		}

		{ //1st list empty
			TR << "Trying with first list empty..." << std::endl;
			utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist2;
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  1, -12.31, 1.23, 8) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  2, -11.41, 2.72, 5) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  3,  -7.93, 3.11, 7) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  4,  -6.33, 4.52, 2) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  5,  -5.99, 1.83, 3) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  6,  -5.54, 0.07, 2) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  7,  -2.12, 1.32, 1) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  8,  -1.43, 0.32, 3) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  9,   1.63, 0.86, 1) ) );
			summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3, 10,   2.12, 1.96, 1) ) );

			utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist1;

			TR << "List 1 before sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
				summarylist1[i]->show(TR); TR.flush();
			}
			TR << std::endl;

			TR << "List 2 before sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist2.size(); i<=imax; ++i ) {
				summarylist2[i]->show(TR); TR.flush();
			}
			TR << std::endl;

			mergesort_jobsummaries_list( summarylist1, summarylist2, SORT_BY_ENERGIES );

			TR << "After sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
				summarylist1[i]->show(TR);  TR.flush();
				if ( i>1 ) TS_ASSERT( summarylist1[i]->pose_energy() >= summarylist1[i-1]->pose_energy() );
			}
			TR << std::endl;
			TS_ASSERT_EQUALS(summarylist1.size(), 10);
		}

		{ //Both lists empty
			TR << "Trying with both lists empty..." << std::endl;
			utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist2;
			utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist1;

			TR << "List 1 before sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
				summarylist1[i]->show(TR); TR.flush();
			}
			TR << std::endl;

			TR << "List 2 before sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist2.size(); i<=imax; ++i ) {
				summarylist2[i]->show(TR); TR.flush();
			}
			TR << std::endl;

			mergesort_jobsummaries_list( summarylist1, summarylist2, SORT_BY_ENERGIES );

			TR << "After sort:" << std::endl;
			for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
				summarylist1[i]->show(TR);  TR.flush();
				if ( i>1 ) TS_ASSERT( summarylist1[i]->pose_energy() >= summarylist1[i-1]->pose_energy() );
			}
			TR << std::endl;
			TS_ASSERT_EQUALS(summarylist1.size(), 0);
		}

	}


	/// @brief Test the sort algorithm used by the MPI version to sort job results by rmsd.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_mergesort_rmsd() {
		TR << "Running SimpleCycpepPredictApplicationTests::test_mergesort_rmsd()." << std::endl;
		TR << "This tests the merge sort algorithm used by the MPI version of the app.  This one tests sorting by RMSD." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;

		using namespace protocols::cyclic_peptide_predict;

		utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist1;
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  1,  -2.12, 0.23, 1) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  2, -11.41, 0.72, 4) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  3,   7.93, 0.83, 1) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  4,  -2.33, 1.52, 2) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  5, -12.99, 1.83, 5) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  6,  -5.54, 2.07, 1) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  7, -12.31, 2.32, 5) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  8,  -1.43, 2.35, 1) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3,  9,   3.63, 3.86, 3) ) );
		summarylist1.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(3, 10,  -5.12, 4.96, 2) ) );

		utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > summarylist2;
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  1, -14.74, 0.21, 6) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  2,  -3.86, 0.78, 2) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  3, -13.32, 0.91, 6) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  4,  -0.28, 1.42, 1) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  5, -19.97, 1.90, 5) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  6,  -8.11, 2.17, 1) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  7, -12.10, 2.33, 5) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  8,  -3.09, 2.37, 3) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5,  9,  -0.06, 2.86, 1) ) );
		summarylist2.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP( new SimpleCycpepPredictApplication_MPI_JobResultsSummary(5, 10,   2.88, 3.46, 2) ) );


		TR << "List 1 before sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
			summarylist1[i]->show(TR); TR.flush();
		}
		TR << std::endl;

		TR << "List 2 before sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist2.size(); i<=imax; ++i ) {
			summarylist2[i]->show(TR); TR.flush();
		}
		TR << std::endl;

		mergesort_jobsummaries_list( summarylist1, summarylist2, SORT_BY_RMSD );

		TR << "After sort:" << std::endl;
		for ( core::Size i=1, imax=summarylist1.size(); i<=imax; ++i ) {
			summarylist1[i]->show(TR);  TR.flush();
			if ( i>1 ) TS_ASSERT( summarylist1[i]->rmsd() >= summarylist1[i-1]->rmsd() );
		}
		TR << std::endl;

	}

	/// @brief Test the read-in of files defining what residues we can and can't design with.
	/// @details This is not scorefunction-dependent, so I am not adding this to the beta_nov15 tests.
	void test_designfile_read() {
		TR << "Running SimpleCycpepPredictApplicationTests::test_designfile_read()." << std::endl;
		TR << "This tests the read-in of files specifying what residues can be used at what positions." << std::endl;
		TR << "For questions, contact Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;

		std::map < core::Size, utility::vector1 < std::string > > allowed_canonicals;
		std::map < core::Size, utility::vector1 < std::string > > allowed_noncanonicals;

		protocols::cyclic_peptide_predict::read_peptide_design_file( "protocols/cyclic_peptide_predict/design_setup_test.txt", allowed_canonicals, allowed_noncanonicals );

		//Check that we've got entries for the correct things:
		for ( core::Size i=0, imax=25; i<=imax; ++i ) {
			if ( i == 0 || i == 1 || i == 3 ) {
				TS_ASSERT_EQUALS( allowed_canonicals.count(i), 1 );
			} else {
				TS_ASSERT_EQUALS( allowed_canonicals.count(i), 0 );
			}
		}

		//Check that we've got the correct number of entries:
		TS_ASSERT_EQUALS( allowed_canonicals.at(0).size(), 2 );
		TS_ASSERT_EQUALS( allowed_canonicals.at(1).size(), 5 );
		TS_ASSERT_EQUALS( allowed_canonicals.at(3).size(), 3 );
		TS_ASSERT_EQUALS( allowed_noncanonicals.at(0).size(), 1 );
		TS_ASSERT_EQUALS( allowed_noncanonicals.at(1).size(), 1 );
		TS_ASSERT_EQUALS( allowed_noncanonicals.at(3).size(), 3 );

		//Check that we've got the correct entries:
		TS_ASSERT_EQUALS( allowed_canonicals.at(0)[1], "ALA"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(0)[2], "GLY"  );
		TS_ASSERT_EQUALS( allowed_noncanonicals.at(0)[1], "DARG"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(1)[1], "ALA"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(1)[2], "ASP"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(1)[3], "GLU"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(1)[4], "ARG"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(1)[5], "LYS"  );
		TS_ASSERT_EQUALS( allowed_noncanonicals.at(1)[1], "NORLEU"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(3)[1], "PHE"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(3)[2], "VAL"  );
		TS_ASSERT_EQUALS( allowed_canonicals.at(3)[3], "LEU"  );
		TS_ASSERT_EQUALS( allowed_noncanonicals.at(3)[1], "DMET"  );
		TS_ASSERT_EQUALS( allowed_noncanonicals.at(3)[2], "DILE"  );
		TS_ASSERT_EQUALS( allowed_noncanonicals.at(3)[3], "DLEU"  );
	}

}; //class GeneralizedKIC_Tests
