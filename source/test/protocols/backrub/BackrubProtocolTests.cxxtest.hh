// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/backrub/BackrubProtocolTests.cxxtest.hh
/// @brief  Unit tests for the protocol level BackrubProtocol
/// @author Kyle Barlow (kb@kylebarlow.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <test/protocols/init_util.hh>
#include <test/util/rosettascripts.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/backrub/BackrubProtocol.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("BackrubProtocolTests");


class BackrubProtocolTests : public CxxTest::TestSuite {
	//Define Variables
private:

	protocols::backrub::BackrubProtocolOP test_instantiation_;
	basic::datacache::DataMap data_;
	Filters_map filters_;
	Movers_map movers_;

public:

	void setUp(){
		core_init();
		test_instantiation_ = protocols::backrub::BackrubProtocolOP( new protocols::backrub::BackrubProtocol() );
	}

	void tearDown(){

	}

	void setUpScripts() {
		prime_Movers( movers_ ); // Adds "null" to movers map
		prime_Filters( filters_ ); // Adds true_filter, false_filter

		// Create dummy "commandline" score function
		core::scoring::ScoreFunctionOP dummy_commandline_sfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
		data_.add( "scorefxns", "commandline", dummy_commandline_sfxn );
	}


	void test_parse_pivots_from_residue_selector(){
		core::select::residue_selector::ResidueSubset residue_subset;
		residue_subset.push_back( true );
		residue_subset.push_back( false );
		residue_subset.push_back( true );
		residue_subset.push_back( true );
		residue_subset.push_back( false );
		test_instantiation_->set_pivots_from_residue_subset( residue_subset );

		utility::vector1<core::Size> pivot_residues = test_instantiation_->get_pivot_residues();

		TS_ASSERT( pivot_residues.size() == 3 );
		TS_ASSERT( pivot_residues[1] == 1 );
		TS_ASSERT( pivot_residues[2] == 3 );
		TS_ASSERT( pivot_residues[3] == 4 );

	}



};
