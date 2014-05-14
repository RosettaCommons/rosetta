// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/prob_util.cxxtest.hh
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <cxxtest/TestSuite.h>

// C/C++ headers
#include <iterator>
#include <vector>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <numeric/prob_util.hh>

#define MARGIN 0.0001

using namespace std;

class ProbUtil : public CxxTest::TestSuite {
 public:
  /// @brief Probabilities in a vector
  vector<double> probs;

  void setUp() {
    const unsigned num_probs = 4;
    const double data[] = { 0.05, 0.15, 0.35, 0.45 };

    vector<double>::iterator i = probs.begin();
    probs.insert(i, data, data + num_probs);
  }

  void tearDown() {
    probs.clear();
  }

  void test_sum() {
    TS_ASSERT_DELTA(1.0, numeric::sum(probs.begin(), probs.end()), MARGIN);
  }

  void test_prenormalized_input() {
    vector<double> tmp(probs.begin(), probs.end());
    numeric::normalize(tmp.begin(), tmp.end());

    // input was already normalized, no change expected
    for (unsigned i = 0; i < probs.size(); ++i)
      TS_ASSERT_DELTA(probs[i], tmp[i], MARGIN);
  }

  void test_normalize() {
    vector<double> x;
    x.push_back(1);
    x.push_back(2);
    x.push_back(3);
    x.push_back(4);

    numeric::normalize(x.begin(), x.end());

    TS_ASSERT_DELTA(0.1, x[0], MARGIN);
    TS_ASSERT_DELTA(0.2, x[1], MARGIN);
    TS_ASSERT_DELTA(0.3, x[2], MARGIN);
    TS_ASSERT_DELTA(0.4, x[3], MARGIN);
  }

  void test_cumulative() {
    numeric::cumulative(probs.begin(), probs.end());

    TS_ASSERT_DELTA(0.05, probs[0], MARGIN);
    TS_ASSERT_DELTA(0.20, probs[1], MARGIN);
    TS_ASSERT_DELTA(0.55, probs[2], MARGIN);
    TS_ASSERT_DELTA(1.00, probs[3], MARGIN);
  }

  void test_product() {
    vector<double> p1;
    p1.push_back(0.5);
    p1.push_back(0.25);
    p1.push_back(0.25);

    vector<double> p2;
    p2.push_back(0.75);
    p2.push_back(0.05);
    p2.push_back(0.20);

    numeric::product(p1.begin(), p1.end(), p2.begin(), p2.end());
    TS_ASSERT_DELTA(0.85714, p1[0], MARGIN);
    TS_ASSERT_DELTA(0.02857, p1[1], MARGIN);
    TS_ASSERT_DELTA(0.11428, p1[2], MARGIN);
  }

  void test_read_probabilities() {
    double data[] = { 0.2, 0.1, 0.4, 0.25, 0.05 };
    utility::vector1<double> expected(data, data + 5);

    utility::vector1<double> observed;
    numeric::read_probabilities_or_die("numeric/sampling.prob", &observed);

		TS_ASSERT_EQUALS( observed.size(), expected.size() );

    for (unsigned i = 1; i <= observed.size(); ++i)
      TS_ASSERT_DELTA(expected[i], observed[i], MARGIN);
  }
};
