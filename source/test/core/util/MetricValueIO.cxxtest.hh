// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/MetricValueIO.cxxtest.hh
/// @brief  test suite for basic::MetricValueIO.cc
/// @author Colin A. Smith


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <basic/MetricValueIO.hh>
#include <basic/MetricValue.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <string>

// C++ headers, for debugging your tests
// AUTO-REMOVED #include <iostream>

//Auto Headers


using std::string;
using core::Real;
using core::Size;
using basic::MetricValueBase;
using basic::MetricValueBaseOP;
using basic::MetricValue;
using basic::check_cast;
using utility::vector1;
using utility::tools::make_vector1;

using basic::write_metric_value;
using basic::read_metric_value;


// --------------- Test Class --------------- //

class MetricValueIOTest : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.

	string real_s_;
	MetricValue<Real> real_mv_;

	string vector1_real_s_;
	MetricValue<vector1<Real> > vector1_real_mv_;

	string vector1_real_empty_s_;
	MetricValue<vector1<Real> > vector1_real_empty_mv_;

	string int_s_;
	MetricValue<int> int_mv_;

	string vector1_int_s_;
	MetricValue<vector1<int> > vector1_int_mv_;

	string size_s_;
	MetricValue<Size> size_mv_;

	string vector1_size_s_;
	MetricValue<vector1<Size> > vector1_size_mv_;

	string bool_s_;
	MetricValue<bool> bool_mv_;

	string vector1_bool_s_;
	MetricValue<vector1<bool> > vector1_bool_mv_;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		real_s_ = "Real 2.56e+10";
		real_mv_.set(2.56e+10);

		vector1_real_s_ = "Real[ 2.3 5.6 1e-10 ]";
		vector1_real_mv_.set(make_vector1<Real>(2.3, 5.6, 1e-10));

		vector1_real_empty_s_ = "Real[ ]";

		int_s_ = "Int -56798";
		int_mv_.set(-56798);

		vector1_int_s_ = "Int[ 4 -89 0 ]";
		vector1_int_mv_.set(make_vector1<int>(4, -89, 0));

		size_s_ = "Size 400";
		size_mv_.set(400);

		vector1_size_s_ = "Size[ 8 10000000 0 ]";
		vector1_size_mv_.set(make_vector1<Size>(8, 10000000, 0));

		bool_s_ = "Bool 1";
		bool_mv_.set(true);

		vector1_bool_s_ = "Bool[ 1 0 1 ]";
		vector1_bool_mv_.set(make_vector1<bool>(true, false, true));
	}

	// Shared finalization goes here.
	void tearDown() {
		// Being a smart pointer, g should be destructed and "free'd" correctly, but a g->delete_everything()
		// could be placed here, if desired.
	}


	// --------------- Test Cases --------------- //

	void test_real() {

		// writing from a MetricValue
		std::ostringstream real_oss1;
		TS_ASSERT(write_metric_value(real_oss1, real_mv_));
		TS_ASSERT_EQUALS(real_s_, real_oss1.str());

		// reading into a MetricValue<Real>
		std::istringstream real_iss1(real_s_);
		MetricValue<Real> real_mv_from_iss1;
		TS_ASSERT(read_metric_value(real_iss1, real_mv_from_iss1));
		TS_ASSERT_EQUALS(real_mv_.value(), real_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream real_iss2(real_s_);
		MetricValueBaseOP real_mvop_from_iss(read_metric_value(real_iss2));
		TS_ASSERT(real_mvop_from_iss);
		TS_ASSERT(check_cast<Real>(real_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			real_mv_.value(),
			static_cast<MetricValue<Real> *> (real_mvop_from_iss.get())->value()
		);

		// reading a Real vector string into a MetricValue<Real> (should fail)
		std::istringstream vector1_real_iss(vector1_real_s_);
		MetricValue<Real> real_mv_from_iss2;
		TS_ASSERT(!read_metric_value(vector1_real_iss, real_mv_from_iss2));
	}

	void test_vector1_real() {

		// writing from a MetricValue
		std::ostringstream vector1_real_oss1;
		TS_ASSERT(write_metric_value(vector1_real_oss1, vector1_real_mv_));
		TS_ASSERT_EQUALS(vector1_real_s_, vector1_real_oss1.str());

		// reading into a MetricValue<vector1<Real> >
		std::istringstream vector1_real_iss1(vector1_real_s_);
		MetricValue<vector1<Real> > vector1_real_mv_from_iss1;
		TS_ASSERT(read_metric_value(vector1_real_iss1, vector1_real_mv_from_iss1));
		TS_ASSERT_EQUALS(vector1_real_mv_.value(), vector1_real_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream vector1_real_iss2(vector1_real_s_);
		MetricValueBaseOP vector1_real_mvop_from_iss(read_metric_value(vector1_real_iss2));
		TS_ASSERT(vector1_real_mvop_from_iss);
		TS_ASSERT(check_cast<vector1<Real> >(vector1_real_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			vector1_real_mv_.value(),
			static_cast<MetricValue<vector1<Real> > *> (vector1_real_mvop_from_iss.get())->value()
		);

		// reading a Real string into a MetricValue<vector1<Real> > (should fail)
		std::istringstream real_iss(real_s_);
		MetricValue<vector1<Real> > vector1_real_mv_from_iss2;
		TS_ASSERT(!read_metric_value(real_iss, vector1_real_mv_from_iss2));
	}

	void test_vector1_real_empty() {

		// writing from a MetricValue
		std::ostringstream vector1_real_empty_oss1;
		TS_ASSERT(write_metric_value(vector1_real_empty_oss1, vector1_real_empty_mv_));
		TS_ASSERT_EQUALS(vector1_real_empty_s_, vector1_real_empty_oss1.str());

		// reading into a MetricValue<vector1<Real> >
		std::istringstream vector1_real_empty_iss1(vector1_real_empty_s_);
		MetricValue<vector1<Real> > vector1_real_empty_mv_from_iss1;
		TS_ASSERT(read_metric_value(vector1_real_empty_iss1, vector1_real_empty_mv_from_iss1));
		TS_ASSERT_EQUALS(vector1_real_empty_mv_.value(), vector1_real_empty_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream vector1_real_empty_iss2(vector1_real_empty_s_);
		MetricValueBaseOP vector1_real_empty_mvop_from_iss(read_metric_value(vector1_real_empty_iss2));
		TS_ASSERT(vector1_real_empty_mvop_from_iss);
		TS_ASSERT(check_cast<vector1<Real> >(vector1_real_empty_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			vector1_real_empty_mv_.value(),
			static_cast<MetricValue<vector1<Real> > *> (vector1_real_empty_mvop_from_iss.get())->value()
		);
	}

	void test_int() {

		// writing from a MetricValue
		std::ostringstream int_oss1;
		TS_ASSERT(write_metric_value(int_oss1, int_mv_));
		TS_ASSERT_EQUALS(int_s_, int_oss1.str());

		// reading into a MetricValue<int>
		std::istringstream int_iss1(int_s_);
		MetricValue<int> int_mv_from_iss1;
		TS_ASSERT(read_metric_value(int_iss1, int_mv_from_iss1));
		TS_ASSERT_EQUALS(int_mv_.value(), int_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream int_iss2(int_s_);
		MetricValueBaseOP int_mvop_from_iss(read_metric_value(int_iss2));
		TS_ASSERT(int_mvop_from_iss);
		TS_ASSERT(check_cast<int>(int_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			int_mv_.value(),
			static_cast<MetricValue<int> *> (int_mvop_from_iss.get())->value()
		);

		// reading a int vector string into a MetricValue<int> (should fail)
		std::istringstream vector1_int_iss(vector1_int_s_);
		MetricValue<int> int_mv_from_iss2;
		TS_ASSERT(!read_metric_value(vector1_int_iss, int_mv_from_iss2));
	}

	void test_vector1_int() {

		// writing from a MetricValue
		std::ostringstream vector1_int_oss1;
		TS_ASSERT(write_metric_value(vector1_int_oss1, vector1_int_mv_));
		TS_ASSERT_EQUALS(vector1_int_s_, vector1_int_oss1.str());

		// reading into a MetricValue<vector1<int> >
		std::istringstream vector1_int_iss1(vector1_int_s_);
		MetricValue<vector1<int> > vector1_int_mv_from_iss1;
		TS_ASSERT(read_metric_value(vector1_int_iss1, vector1_int_mv_from_iss1));
		TS_ASSERT_EQUALS(vector1_int_mv_.value(), vector1_int_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream vector1_int_iss2(vector1_int_s_);
		MetricValueBaseOP vector1_int_mvop_from_iss(read_metric_value(vector1_int_iss2));
		TS_ASSERT(vector1_int_mvop_from_iss);
		TS_ASSERT(check_cast<vector1<int> >(vector1_int_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			vector1_int_mv_.value(),
			static_cast<MetricValue<vector1<int> > *> (vector1_int_mvop_from_iss.get())->value()
		);

		// reading a int string into a MetricValue<vector1<int> > (should fail)
		std::istringstream int_iss(int_s_);
		MetricValue<vector1<int> > vector1_int_mv_from_iss2;
		TS_ASSERT(!read_metric_value(int_iss, vector1_int_mv_from_iss2));
	}

	void test_size() {

		// writing from a MetricValue
		std::ostringstream size_oss1;
		TS_ASSERT(write_metric_value(size_oss1, size_mv_));
		TS_ASSERT_EQUALS(size_s_, size_oss1.str());

		// reading into a MetricValue<Size>
		std::istringstream size_iss1(size_s_);
		MetricValue<Size> size_mv_from_iss1;
		TS_ASSERT(read_metric_value(size_iss1, size_mv_from_iss1));
		TS_ASSERT_EQUALS(size_mv_.value(), size_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream size_iss2(size_s_);
		MetricValueBaseOP size_mvop_from_iss(read_metric_value(size_iss2));
		TS_ASSERT(size_mvop_from_iss);
		TS_ASSERT(check_cast<Size>(size_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			size_mv_.value(),
			static_cast<MetricValue<Size> *> (size_mvop_from_iss.get())->value()
		);

		// reading a Size vector string into a MetricValue<Size> (should fail)
		std::istringstream vector1_size_iss(vector1_size_s_);
		MetricValue<Size> size_mv_from_iss2;
		TS_ASSERT(!read_metric_value(vector1_size_iss, size_mv_from_iss2));
	}

	void test_vector1_size() {

		// writing from a MetricValue
		std::ostringstream vector1_size_oss1;
		TS_ASSERT(write_metric_value(vector1_size_oss1, vector1_size_mv_));
		TS_ASSERT_EQUALS(vector1_size_s_, vector1_size_oss1.str());

		// reading into a MetricValue<vector1<Size> >
		std::istringstream vector1_size_iss1(vector1_size_s_);
		MetricValue<vector1<Size> > vector1_size_mv_from_iss1;
		TS_ASSERT(read_metric_value(vector1_size_iss1, vector1_size_mv_from_iss1));
		TS_ASSERT_EQUALS(vector1_size_mv_.value(), vector1_size_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream vector1_size_iss2(vector1_size_s_);
		MetricValueBaseOP vector1_size_mvop_from_iss(read_metric_value(vector1_size_iss2));
		TS_ASSERT(vector1_size_mvop_from_iss);
		TS_ASSERT(check_cast<vector1<Size> >(vector1_size_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			vector1_size_mv_.value(),
			static_cast<MetricValue<vector1<Size> > *> (vector1_size_mvop_from_iss.get())->value()
		);

		// reading a Size string into a MetricValue<vector1<Size> > (should fail)
		std::istringstream size_iss(size_s_);
		MetricValue<vector1<Size> > vector1_size_mv_from_iss2;
		TS_ASSERT(!read_metric_value(size_iss, vector1_size_mv_from_iss2));
	}

	void test_bool() {

		// writing from a MetricValue
		std::ostringstream bool_oss1;
		TS_ASSERT(write_metric_value(bool_oss1, bool_mv_));
		TS_ASSERT_EQUALS(bool_s_, bool_oss1.str());

		// reading into a MetricValue<bool>
		std::istringstream bool_iss1(bool_s_);
		MetricValue<bool> bool_mv_from_iss1;
		TS_ASSERT(read_metric_value(bool_iss1, bool_mv_from_iss1));
		TS_ASSERT_EQUALS(bool_mv_.value(), bool_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream bool_iss2(bool_s_);
		MetricValueBaseOP bool_mvop_from_iss(read_metric_value(bool_iss2));
		TS_ASSERT(bool_mvop_from_iss);
		TS_ASSERT(check_cast<bool>(bool_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			bool_mv_.value(),
			static_cast<MetricValue<bool> *> (bool_mvop_from_iss.get())->value()
		);

		// reading a bool vector string into a MetricValue<bool> (should fail)
		std::istringstream vector1_bool_iss(vector1_bool_s_);
		MetricValue<bool> bool_mv_from_iss2;
		TS_ASSERT(!read_metric_value(vector1_bool_iss, bool_mv_from_iss2));
	}

	void test_vector1_bool() {

		// writing from a MetricValue
		std::ostringstream vector1_bool_oss1;
		TS_ASSERT(write_metric_value(vector1_bool_oss1, vector1_bool_mv_));
		TS_ASSERT_EQUALS(vector1_bool_s_, vector1_bool_oss1.str());

		// reading into a MetricValue<vector1<bool> >
		std::istringstream vector1_bool_iss1(vector1_bool_s_);
		MetricValue<vector1<bool> > vector1_bool_mv_from_iss1;
		TS_ASSERT(read_metric_value(vector1_bool_iss1, vector1_bool_mv_from_iss1));
		TS_ASSERT_EQUALS(vector1_bool_mv_.value(), vector1_bool_mv_from_iss1.value());

		// reading into a MetricValueBaseOP
		std::istringstream vector1_bool_iss2(vector1_bool_s_);
		MetricValueBaseOP vector1_bool_mvop_from_iss(read_metric_value(vector1_bool_iss2));
		TS_ASSERT(vector1_bool_mvop_from_iss);
		TS_ASSERT(check_cast<vector1<bool> >(vector1_bool_mvop_from_iss.get()));
		TS_ASSERT_EQUALS(
			vector1_bool_mv_.value(),
			static_cast<MetricValue<vector1<bool> > *> (vector1_bool_mvop_from_iss.get())->value()
		);

		// reading a bool string into a MetricValue<vector1<bool> > (should fail)
		std::istringstream bool_iss(bool_s_);
		MetricValue<vector1<bool> > vector1_bool_mv_from_iss2;
		TS_ASSERT(!read_metric_value(bool_iss, vector1_bool_mv_from_iss2));
	}

	void test_unsupported_value() {

		std::ostringstream oss;
		MetricValue<std::string> string_mv;
		TS_ASSERT(!write_metric_value(oss, string_mv));
		MetricValue<vector1<std::string> > vector1_string_mv;
		TS_ASSERT(!write_metric_value(oss, vector1_string_mv));
	}

	void test_bad_formatting() {

		// all of these should fail
		std::istringstream iss1("Bool ");
		TS_ASSERT(!read_metric_value(iss1));
		std::istringstream iss2("Bool true"); // bools are written as 1 and 0!
		TS_ASSERT(!read_metric_value(iss2));
		std::istringstream iss3("Real[ 1 2 5");
		TS_ASSERT(!read_metric_value(iss3));
		std::istringstream iss4("Size[0 34 5 ]");
		TS_ASSERT(!read_metric_value(iss4));
		std::istringstream iss5("Size[ 0 34 5]");
		TS_ASSERT(!read_metric_value(iss5));
		std::istringstream iss6("Size [ 0 34 5 ]");
		TS_ASSERT(!read_metric_value(iss6));
		std::istringstream iss7("String \"Foo\""); // strings aren't supported, this test should be updated when they are
		TS_ASSERT(!read_metric_value(iss7));
		std::istringstream iss8("String[ \"Foo\" \"Bar\" ]");
		TS_ASSERT(!read_metric_value(iss8));
		std::istringstream iss9("");
		TS_ASSERT(!read_metric_value(iss9));
		std::istringstream iss10(" ");
		TS_ASSERT(!read_metric_value(iss10));
		std::istringstream iss11("45");
		TS_ASSERT(!read_metric_value(iss11));
	}

};


