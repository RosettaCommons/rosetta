// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/bin_transitions/BinTransitionCalculator.cxxtest.hh
/// @brief  Unit tests for the BinTransitionCalculator class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers:
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// BinTransitionCalculator headers:
#include <core/scoring/bin_transitions/BinTransitionCalculator.fwd.hh>
#include <core/scoring/bin_transitions/BinTransitionCalculator.hh>
#include <core/scoring/bin_transitions/BinTransitionData.fwd.hh>
#include <core/scoring/bin_transitions/BinTransitionData.hh>

// Other Rosetta libraries:
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

// C++ headers
#include <utility>
#include <iostream>
#include <iomanip>

// --------------- Test Class --------------- //

static thread_local basic::Tracer TR("core.scoring.bin_transitions.BinTransitionCalculator.cxxtest");

class BinTransitionCalculatorTests : public CxxTest::TestSuite {

private:

	core::pose::PoseOP pose_;

public:

	void setUp() {
		core_init();
		pose_ = core::pose::PoseOP( new core::pose::Pose );
		core::pose::make_pose_from_sequence(*pose_, "AAAAGAPAGPAA", "fa_standard", false);
	}

	void tearDown() {
	}

	/// @brief Test the BinTransitionCalculator's read from file.
	/// @details This reads the ABBA.bin_params file and checks that
	/// certain things were set properly.
	void test_BinTransitionCalculator_read()
	{
		using namespace core::scoring::bin_transitions;

		BinTransitionCalculatorOP bt( new BinTransitionCalculator );

		bt->load_bin_params( "ABBA" );
		std::string summary(bt->summarize_stored_data(false)); //Get a summary.
		if ( TR.visible() ) {
			TR << summary << std::endl; //Write out the summary.
			TR.flush();
		}
		return;
	}

	/// @brief Test the BinTransitionData class trim_subbin_edges_and_rescale_subbin() function.
	///
	void test_BinTransitionData_trim_subbin_edges_and_rescale_subbin()
	{
		using namespace core::scoring::bin_transitions;
		using namespace std;

		BinTransitionDataOP btd( new BinTransitionData);
		std::ostringstream outstream;

		core::Real curbin_val(2.0);
		std::pair <core::Real, core::Real> tors_range = make_pair(35, 45);
		core::Real phimin(30), phimax(40);
		outstream << std::endl;
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		btd->trim_subbin_edges_and_rescale_subbin( curbin_val, tors_range, phimin, phimax);
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		outstream << std::endl;
		TS_ASSERT_DELTA(curbin_val, 1.0, 0.00001  );
		TS_ASSERT_DELTA(tors_range.first, 35, 0.00001);
		TS_ASSERT_DELTA(tors_range.second, 40, 0.00001);

		curbin_val=2.0;
		tors_range=make_pair(35,45);
		phimin=40; phimax=50;
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		btd->trim_subbin_edges_and_rescale_subbin( curbin_val, tors_range, phimin, phimax);
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		outstream << std::endl;
		TS_ASSERT_DELTA(curbin_val, 1.0, 0.00001  );
		TS_ASSERT_DELTA(tors_range.first, 40, 0.00001);
		TS_ASSERT_DELTA(tors_range.second, 45, 0.00001);

		curbin_val=2.0;
		tors_range=make_pair(175,-175);
		phimin=170; phimax=180;
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		btd->trim_subbin_edges_and_rescale_subbin( curbin_val, tors_range, phimin, phimax);
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		outstream << std::endl;
		TS_ASSERT_DELTA(curbin_val, 1.0, 0.00001  );
		TS_ASSERT_DELTA(tors_range.first, 175, 0.00001);
		TS_ASSERT_DELTA(tors_range.second, 180, 0.00001);

		curbin_val=2.0;
		tors_range=make_pair(175,-175);
		phimin=172; phimax=-178;
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		btd->trim_subbin_edges_and_rescale_subbin( curbin_val, tors_range, phimin, phimax);
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		outstream << std::endl;
		TS_ASSERT_DELTA(curbin_val, 1.4, 0.00001  );
		TS_ASSERT_DELTA(tors_range.first, 175, 0.00001);
		TS_ASSERT_DELTA(tors_range.second, -178, 0.00001);

		curbin_val=2.0;
		tors_range=make_pair(171,-179);
		phimin=175; phimax=-175;
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		btd->trim_subbin_edges_and_rescale_subbin( curbin_val, tors_range, phimin, phimax);
		outstream << "curbin_val=" << curbin_val << "\ttors_range=" << tors_range.first << "," << tors_range.second << "\tphimin=" << phimin << "\tphimax=" << phimax << std::endl;
		outstream << std::endl;
		TS_ASSERT_DELTA(curbin_val, 1.2, 0.00001  );
		TS_ASSERT_DELTA(tors_range.first, 175, 0.00001);
		TS_ASSERT_DELTA(tors_range.second, -179, 0.00001);

		TR << outstream.str();
		TR.flush();

		return;
	}

	/// @brief Test the probability transition calculator.
	/// @details This test will fail if the bin transition probability data change.
	void test_BinTransitionCalculator_p_iplus1_given_i()
	{ //AAAAGAPAGPAA
		using namespace core::scoring::bin_transitions;
		using namespace std;

		BinTransitionCalculatorOP bt( new BinTransitionCalculator );
		bt->load_bin_params( "ABBA" );

		//A->A
		core::pose::PoseOP temppose = pose_->clone();
		for ( core::Size i=1, imax=temppose->n_residue(); i<=imax; ++i ) {
			if ( i>1 ) temppose->set_phi( i, -64.8 );
			if ( i<imax ) {
				temppose->set_psi( i, -41.0 );
				temppose->set_omega( i, 180.0 );
			}
		}

		core::Real prob1(0);
		core::Real prob2(0);
		core::Real prob3(0);
		bool calc1_success(false);
		bool calc2_success(false);
		bool calc3_success(false);
		core::Real prob1_expected(0.8416743620);
		core::Real prob2_expected(0.5110111721);
		core::Real prob3_expected(0.9446337308);

		calc1_success = bt->p_iplus1_given_i( temppose->residue(2), temppose->residue(3), prob1 ); //ala->ala, A->A
		calc2_success = bt->p_iplus1_given_i( temppose->residue(5), temppose->residue(6), prob2 ); //gly->ala, A->A
		calc3_success = bt->p_i_given_iplus1( temppose->residue(9), temppose->residue(10), prob3 ); //gly->pro, A->A

		//B->Aprime
		for ( core::Size i=1, imax=temppose->n_residue(); i<=imax; ++i ) {
			if ( i>1 ) temppose->set_phi( i, -135 );
			if ( i<imax ) {
				temppose->set_psi( i, 135.0 );
				temppose->set_omega( i, 180.0 );
			}
		}
		temppose->set_phi(3, 64.8);
		temppose->set_psi(3, 41.0);
		temppose->set_phi(6, 64.8);
		temppose->set_psi(6, 41.0);
		temppose->set_phi(10, 64.8);
		temppose->set_psi(10, 41.0);

		core::Real prob4(1);
		core::Real prob5(1);
		core::Real prob6(1);
		bool calc4_success(false);
		bool calc5_success(false);
		bool calc6_success(false);
		core::Real prob4_expected(0.0330397173);
		core::Real prob5_expected(0.4578725398);
		core::Real prob6_expected(0.6347031963);

		calc4_success = bt->p_iplus1_given_i( temppose->residue(2), temppose->residue(3), prob4 ); //ala->ala, B->Aprime
		calc5_success = bt->p_iplus1_given_i( temppose->residue(5), temppose->residue(6), prob5 ); //gly->ala, B->Aprime
		calc6_success = bt->p_i_given_iplus1( temppose->residue(9), temppose->residue(10), prob6 ); //gly->pro, B->Aprime

		if ( TR.visible() ) {
			TR << "Residues\tA-A\tB-A'" << std::endl;
			TR << "ala->ala\t" << setprecision(6) << prob1 << "\t" << setprecision(6) << prob4 << std::endl;
			TR << "gly->ala\t" << setprecision(6) << prob2 << "\t" << setprecision(6) << prob5 << std::endl;
			TR << "gly<-pro\t" << setprecision(6) << prob3 << "\t" << setprecision(6) << prob6 << std::endl;
		}

		//Confirm that bin transition probabilities were found for all of these:
		TS_ASSERT(calc1_success);
		TS_ASSERT(calc2_success);
		TS_ASSERT(calc3_success);
		TS_ASSERT(calc4_success);
		TS_ASSERT(calc5_success);
		TS_ASSERT(calc6_success);

		//Confirm that the values are right. This will have to be updated if the database file changes.
		TS_ASSERT_DELTA( prob1, prob1_expected, 1e-8 );
		TS_ASSERT_DELTA( prob2, prob2_expected, 1e-8 );
		TS_ASSERT_DELTA( prob3, prob3_expected, 1e-8 );
		TS_ASSERT_DELTA( prob4, prob4_expected, 1e-8 );
		TS_ASSERT_DELTA( prob5, prob5_expected, 1e-8 );
		TS_ASSERT_DELTA( prob6, prob6_expected, 1e-8 );

		return;
	}

}; //class BinTransitionCalculatorTests
