// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Field.cxxtest.hh
/// @brief  test suite for classes associated with core::io::pdb::Field
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/pdb/Field.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/pdb/build_pose_as_is.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>


// C++ headers
#include <sstream>
#include <vector>

static basic::Tracer TR("core.io.pdb.Field.cxxtest");

using namespace core;

class FieldTests : public CxxTest::TestSuite
{

public:
	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	void test_field_basic() {
		using core::io::pdb::Field;
		using core::Size;
		using std::string;
		using std::stringstream;

		Size start_first( 1 ), end_first( 5 );
		Size start_second( 6 ), end_second( 11 );
		Size start_third( 12 ), end_third( 16 );
		string record_string( "FIRSTSECONDTHIRD" );

		Field first( start_first, end_first );
		Field second( start_second, end_second );
		Field third( start_third, end_third );

		first.set_value_from_pdb_line( record_string );
		second.set_value_from_pdb_line( record_string );
		third.set_value_from_pdb_line( record_string );

		TS_ASSERT_EQUALS( first.value, "FIRST" );
		TS_ASSERT_EQUALS( second.value, "SECOND" );
		TS_ASSERT_EQUALS( third.value, "THIRD" );

		stringstream first_out, second_out, third_out;
		first_out << first;
		second_out << second;
		third_out << third;

		TS_ASSERT_EQUALS( first_out.str(), "[1, 5]=FIRST" );
		TS_ASSERT_EQUALS( second_out.str(), "[6, 11]=SECOND" );
		TS_ASSERT_EQUALS( third_out.str(), "[12, 16]=THIRD" );
	}
};

