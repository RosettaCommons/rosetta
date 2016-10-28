// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/enzymes/consensus_sequence_parsers.cxxtest.hh
/// @brief   Test suite for consensus sequence parsersing helper functions.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/enzymes/consensus_sequence_parsers.hh>

// Utility header
#include <utility/vector1.hh>

// C++ header
#include <string>


class ConsensusSequenceParserTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that peptide consensus sequences are parsed properly.
	void test_get_3_letter_codes_from_peptide_consensus_sequence()
	{
		using namespace std;
		using namespace utility;

		TS_TRACE( "Testing get_3_letter_codes_from_peptide_consensus_sequence()..." );

		string const sequence( "HE(S/H/E)ISX" );
		vector1< vector1< string > > consensus_residues;

		consensus_residues = core::enzymes::get_3_letter_codes_from_peptide_consensus_sequence( sequence );

		TS_ASSERT_EQUALS( consensus_residues.size(), 6 );
		TS_ASSERT_EQUALS( consensus_residues[ 1 ].size(), 1 );
		TS_ASSERT_EQUALS( consensus_residues[ 1 ][ 1 ], "HIS" );
		TS_ASSERT_EQUALS( consensus_residues[ 3 ].size(), 3 );
		TS_ASSERT_EQUALS( consensus_residues[ 3 ][ 1 ], "SER" );
		TS_ASSERT_EQUALS( consensus_residues[ 3 ][ 2 ], "HIS" );
		TS_ASSERT_EQUALS( consensus_residues[ 3 ][ 3 ], "GLU" );
		TS_ASSERT_EQUALS( consensus_residues[ 4 ].size(), 1 );
		TS_ASSERT_EQUALS( consensus_residues[ 6 ].size(), 20 );
	}
};  // class ConsensusSequenceParserTests
