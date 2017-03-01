// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/UTools.cxxtest.hh
/// @brief  test suite for test/UTools functions.
/// @author Sergey Lyskov

#include <test/UTools.hh>

// Test headers
#include <cxxtest/TestSuite.h>

class UToolsTests : public CxxTest::TestSuite
{
public:
	UToolsTests() {}

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// ------------------------------------------ //
	/// @brief test isFloatNumber function
	void test_isFloatNumber() {
		TS_ASSERT_EQUALS(true, false);

		double num;
		unsigned int end_pos;
		bool res;

		res = test::utools::isFloatNumber("", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber("", 1, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" ", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" ", 1, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" ", 2, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" . ", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" +- ", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" 1", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, 1);
		TS_ASSERT_EQUALS(end_pos, 2);

		res = test::utools::isFloatNumber(" 1 ", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, 1);
		TS_ASSERT_EQUALS(end_pos, 2);

		res = test::utools::isFloatNumber("  1", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" +.", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" +.e0", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber("+.0e0", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" +.0e0", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, 0);
		TS_ASSERT_EQUALS(end_pos, 6);

		res = test::utools::isFloatNumber(" -.1e 111", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" -.1e0 111", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, -.1);
		TS_ASSERT_EQUALS(end_pos, 6);

		res = test::utools::isFloatNumber(" 99.23AX", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" 99.23 AX", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, 99.23);
		TS_ASSERT_EQUALS(end_pos, 6);

		res = test::utools::isFloatNumber(" 99.23 AX", 1, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" 99.23 AX", 4, num, end_pos);
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isFloatNumber(" +12.345678e-10", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, +12.345678e-10);
		TS_ASSERT_EQUALS(end_pos, 15);

		res = test::utools::isFloatNumber(" -12.345678E+10 ___", 0, num, end_pos);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, -12.345678e10);
		TS_ASSERT_EQUALS(end_pos, 15);

		//TS_ASSERT_DELTA(x,y,delta).
	}


	// ------------------------------------------ //
	/// @brief test isMarkup function
	void test_isMarkup() {
		double num;
		unsigned int end_pos;
		bool res;

		res = test::utools::isMarkup("something(1.23)", 0, num, end_pos, "MARKUP(");
		TS_ASSERT_EQUALS(res, false);

		res = test::utools::isMarkup(" MARKUP(1.23)", 0, num, end_pos, "MARKUP(");
		TS_ASSERT_EQUALS(res, false);


		res = test::utools::isMarkup("MARKUP(1.23e-9)", 0, num, end_pos, "MARKUP(");
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, 1.23e-9);
		TS_ASSERT_EQUALS(end_pos, 15);

		res = test::utools::isMarkup("MARKUP(1.23e-9)55566 aa ", 0, num, end_pos, "MARKUP(");
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(num, 1.23e-9);
		TS_ASSERT_EQUALS(end_pos, 15);
	}


	// ------------------------------------------ //
	/// @brief test isEq function
	void test_isEq() {
		double abs_tolerance, rel_tolerance;
		std::string error_message;
		bool res;


		abs_tolerance=0;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("", "", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(error_message, "");

		abs_tolerance=0;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("1 ",
								 "1 ", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(error_message, "");

		abs_tolerance=0;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("set_rel_tolerance(.1) 1.091 aa 5",
								 " 1 aa 4.561", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(abs_tolerance, 0);
		TS_ASSERT_EQUALS(rel_tolerance, .1);
		TS_ASSERT_EQUALS(error_message, "");

		abs_tolerance=0;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("1 set_abs_tolerance(.112)",
								 "1 set_abs_tolerance(.112)", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(abs_tolerance, .112);
		TS_ASSERT_EQUALS(rel_tolerance, 0);
		TS_ASSERT_EQUALS(error_message, "");

		abs_tolerance=0;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("1 set_abs_tolerance(.112)set_rel_tolerance(98.765)",
								 "1 set_abs_tolerance(.2)set_rel_tolerance(.3)", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(abs_tolerance, .112);
		TS_ASSERT_EQUALS(rel_tolerance, 98.765);
		TS_ASSERT_EQUALS(error_message, "");


		abs_tolerance=0;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("1.2e-10 set_abs_tolerance(.091) 1 ",
								 "1.2e-10 set_abs_tolerance(.2) 1.09 ", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(abs_tolerance, .091);
		TS_ASSERT_EQUALS(rel_tolerance, 0);
		TS_ASSERT_EQUALS(error_message, "");


		abs_tolerance=0;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("aa b 1 d set_rel_tolerance(.1) 12 ",
								 "aa b 1 d set_rel_tolerance(.2) 13 ", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, true);
		TS_ASSERT_EQUALS(rel_tolerance, .1);
		TS_ASSERT_EQUALS(abs_tolerance, 0);
		TS_ASSERT_EQUALS(error_message, "");


		abs_tolerance=1;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("aa b 1 d set_abs_tolerance(.0) 12 ",
								 "aa b 1 d set_abs_tolerance(.2) 13 ", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, false);
		TS_ASSERT_DIFFERS(error_message, "");

		abs_tolerance=1;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("aa b x 1 d ",
								 "aa b 1 d ", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, false);
		TS_ASSERT_DIFFERS(error_message, "");


		abs_tolerance=.1;  rel_tolerance=0;  error_message="";
		res = test::utools::isEq("aa b 10 d ",
								 "aa b 10.1001 d ", abs_tolerance, rel_tolerance, error_message);
		TS_ASSERT_EQUALS(res, false);
		TS_ASSERT_DIFFERS(error_message, "");
	}
};
