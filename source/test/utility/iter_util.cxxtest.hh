// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/iter_util.cxxtest.hh
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <cxxtest/TestSuite.h>

// Project headers
#include <utility/iter_util.hh>

// C/C++ headers
#include <iterator>
#include <set>

using namespace std;

class IterUtilTest : public CxxTest::TestSuite {
 public:
  set<int> values;

  void setUp() {
    int data[] = { 1, 2, 5 };
    values.insert(data, data + 3);
  }

  void tearDown() {
    values.clear();
  }

  /// @brief Values contained in range
  void test_find_closest_value_included() {
    TS_ASSERT_EQUALS(1, *utility::find_closest(values.begin(), values.end(), 1));
    TS_ASSERT_EQUALS(2, *utility::find_closest(values.begin(), values.end(), 2));
    TS_ASSERT_EQUALS(5, *utility::find_closest(values.begin(), values.end(), 5));
  }

  /// @brief Values not contained in range
  void test_find_closest_value_excluded() {
    TS_ASSERT_EQUALS(2, *utility::find_closest(values.begin(), values.end(), 3));
    TS_ASSERT_EQUALS(5, *utility::find_closest(values.begin(), values.end(), 4));
  }

  /// @brief Value smaller than anything in the range
  void test_find_closest_smallest() {
    TS_ASSERT_EQUALS(1, *utility::find_closest(values.begin(), values.end(), 0));
  }

  /// @brief Value larger than anything in the range
  void test_find_closest_largest() {
    TS_ASSERT_EQUALS(5, *utility::find_closest(values.begin(), values.end(), 7));
  }

  /// @brief Value equidistant from two values in the range
  void test_find_closest_rounding() {
    set<int> values;
    values.insert(1);
    values.insert(3);

    TS_ASSERT_EQUALS(1, *utility::find_closest(values.begin(), values.end(), 2));
  }
};
