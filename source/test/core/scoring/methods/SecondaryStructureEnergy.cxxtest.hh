// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SecondaryStructureEnergy.cxxtest.hh
/// @brief  test suite for centroid energies that depend on secondary structure
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/SecondaryStructureEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <core/io/silent/SilentFileData.hh>

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class SecondaryStructureEnergy_Tests : public CxxTest::TestSuite {

public:

	Pose pose;
	core::io::silent::SilentFileData sfd;
	core::chemical::ResidueTypeSetCAP rsd_set;
	core::scoring::ScoreFunctionOP scorefxn;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		using namespace std;

		//extern int command_line_argc; extern char ** command_line_argv;
		using namespace core;
		using namespace core::scoring;
		core_init();

		// correct answers taken from rosetta++ v19429
		sfd.read_file( "core/scoring/methods/score3_in.silent_out" );
		// ss energy calculations use centroid residue types
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		scorefxn = new core::scoring::ScoreFunction;
		scorefxn->set_weight( hs_pair, 1.0 );
		scorefxn->set_weight( ss_pair, 1.0 );
		scorefxn->set_weight( rsigma,  1.0 );
		scorefxn->set_weight( sheet,   1.0 );
	}

	// Shared finalization goes here.
	void tearDown() {
		scorefxn = 0;
	}

	// --------------- Test Cases --------------- //
	void test_eval_energy() {
		// TOLERATED_ERROR is set to the largest deviation from calculating contact order in mini and r++
		// as of 8/13/08.
		float const TOLERATED_ERROR = 5.0;

		float hs_pair, ss_pair, rsigma, sheet; // scores calculated in mini
		float silent_hs_pair, silent_ss_pair, silent_rsigma, silent_sheet; // scores from r++ silent-file
		EnergyMap emap;
		int ngood = 0;
		for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
			iter->fill_pose( pose, *rsd_set );

			//emap.zero();
			// scores from silent-file
			silent_hs_pair = iter->get_energy( "hs_pair" );
			silent_ss_pair = iter->get_energy( "ss_pair" );
			silent_rsigma  = iter->get_energy( "rsigma"  );
			silent_sheet   = iter->get_energy( "sheet"   );

			// calculated scores in mini
			//ss_energy.setup_for_scoring( pose, *scorefxn );
			//ss_energy.finalize_total_energy( pose, *scorefxn, emap );
			(*scorefxn)(pose);
			emap = pose.energies().total_energies();
			hs_pair = emap[ core::scoring::hs_pair ];
			ss_pair = emap[ core::scoring::ss_pair ];
			rsigma  = emap[ core::scoring::rsigma  ];
			sheet   = emap[ core::scoring::sheet   ];

			TS_ASSERT_DELTA( hs_pair, silent_hs_pair, TOLERATED_ERROR );
			TS_ASSERT_DELTA( ss_pair, silent_ss_pair, TOLERATED_ERROR );
			TS_ASSERT_DELTA( rsigma,  silent_rsigma,  TOLERATED_ERROR );
			TS_ASSERT_DELTA( sheet,   silent_sheet,   TOLERATED_ERROR );
			// emap.show_if_nonzero_weight( std::cout, scorefxn->weights() );
			if ( std::abs( hs_pair - silent_hs_pair ) <= TOLERATED_ERROR ) {
				ngood++;
			}
			//std::cout << "decoy_tag = " << iter->decoy_tag() << std::endl;
		}

		//std::cout << std::endl << "found " << ngood << " hs_pair scores that passed threshold!" << std::endl;
	} // test_eval_energy()
}; // SecondaryStructureEnergy_Tests


