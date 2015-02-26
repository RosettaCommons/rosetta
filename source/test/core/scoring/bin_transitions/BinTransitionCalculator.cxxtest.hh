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

// C++ headers
#include <utility>
#include <iostream>

// --------------- Test Class --------------- //

static thread_local basic::Tracer TR("core.scoring.bin_transitions.BinTransitionCalculator.cxxtest");

class BinTransitionCalculatorTests : public CxxTest::TestSuite {

private:

public:

	void setUp() {
		core_init();

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
		if(TR.visible()) {
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

}; //class BinTransitionCalculator_Tests
