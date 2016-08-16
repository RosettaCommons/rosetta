// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/minmax.cxxtest.hh
/// @brief  test suite for utility::minmax.hh
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// C/C++ headers
#include <vector>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <utility/minmax.hh>

namespace {

using std::vector;
using utility::vector1;

class MinMaxTest : public CxxTest::TestSuite {
 public:
	vector<double> vector_;
	vector1<double> vector1_;

  void setUp() {
    for (int i = 100; i > 0; --i) {
			vector_.push_back(i);
			vector1_.push_back(i);
		}
  }

  void tearDown() {
		vector_.clear();
		vector1_.clear();
	}

	// --vector-- //

	void test_argmin_std_vector() {
		TS_ASSERT_EQUALS(vector_.size() - 1, utility::argmin(vector_));
	}

	void test_argmin_std_vector_empty_input() {
		vector<double> empty;
		TS_ASSERT_EQUALS(-1, utility::argmin(empty));
	}

	void test_argmax_std_vector() {
		TS_ASSERT_EQUALS(0, utility::argmax(vector_));
	}

	void test_argmax_std_vector_empty_input() {
		vector<double> empty;
		TS_ASSERT_EQUALS(-1, utility::argmax(empty));
	}

	// -- vector1 -- //

	void test_argmin_vector1() {
		TS_ASSERT_EQUALS(vector1_.size(), utility::argmin(vector1_));
	}

	void test_argmin_vector1_empty_input() {
		vector1<double> empty;
		TS_ASSERT_EQUALS(0, utility::argmin(empty));
	}

	void test_argmax_vector1() {
		TS_ASSERT_EQUALS(1, utility::argmax(vector1_));
	}

	void test_argmax_vector1_empty_input() {
		vector1<double> empty;
		TS_ASSERT_EQUALS(0, utility::argmax(empty));
	}
};

}  // anonymous namespace
