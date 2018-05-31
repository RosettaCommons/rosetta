// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/enzymes/database_io.cxxtest.hh
/// @brief   Test suite for enzyme data input/output from the database.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/enzymes/database_io.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "EnzymeDatabaseIOTests" );

class EnzymeDatabaseIOTests : public CxxTest::TestSuite {
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
	void test_read_enzyme_data_from_file()
	{
		using namespace core::enzymes;

		TR << "Testing read_enzyme_data_from_file()..." << std::endl;

		EnzymeData const rosettase_data( read_enzyme_data_from_file( "core/enzymes/rosettase" ) );

		TS_ASSERT_EQUALS( rosettase_data.consensus_sequence, "TARGET" );
		TS_ASSERT_EQUALS( rosettase_data.cs_type, AA );
		TS_ASSERT_EQUALS( rosettase_data.cs_resnum, 4 );
		TS_ASSERT_EQUALS( rosettase_data.atom_to_modify, "CA" );
		TS_ASSERT_EQUALS( rosettase_data.efficiency, 1.00 );
		TS_ASSERT_EQUALS( rosettase_data.second_substrates_or_byproducts.size(), 5 );
		TS_ASSERT_EQUALS( rosettase_data.second_substrates_or_byproducts[ 1 ], "ARROW" );
		TS_ASSERT_EQUALS( rosettase_data.second_substrates_or_byproducts[ 3 ], "BOLT" );
		TS_ASSERT_EQUALS( rosettase_data.second_substrates_or_byproducts[ 5 ], "PROTONTORPEDO" );
	}
};  // class DatabaseIOTests
