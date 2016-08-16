// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   test/utility/signals/PausableSignalHub.cxxtest.hh
/// @brief  test for utility::signals::PausableSignalHub
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <utility/signals/PausableSignalHub.hh>

// C++ headers
#include <vector>


// here must define a different namespace than in SignalHub.cxxtest.hh/BufferedSignalHub.cxxtest.hh
// otherwise we run into strange crosstalk situation and tests begin to fail
namespace sig3 {


// Event
struct Event {};

} // namespace sig3


// --------------- Test Class --------------- //

class PausableSignalHubTests : public CxxTest::TestSuite {


public: // setup


  // shared initialization
  void setUp() {
  }


  // shared finalization
  void tearDown() {
  }


  // --------------- Test Cases --------------- //

	/// @brief test observer pausing, we mainly test that this class compiles ok
	void test_PausableSignalHub_pausing() {
		using namespace sig3;

		utility::signals::PausableSignalHub< void, Event > hub;

		// we don't actually pause here, because it would interrupt the
		// running of the unit tests
		hub.pause();
		TS_ASSERT( hub.pausing() );

		hub.unpause();
		TS_ASSERT( !hub.pausing() );
	}


}; // class BufferedSignalHubTests

