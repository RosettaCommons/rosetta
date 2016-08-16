// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/random/WeightedReservoirSampler.cxxtest.hh
/// @brief test suite for numeric/random/WeightedReservoirSampler.hh
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <numeric/random/random.hh>
#include <numeric/random/WeightedReservoirSampler.hh>

namespace {

using numeric::random::WeightedReservoirSampler;
using utility::vector1;

class WeightedReservoirSamplerTest : public CxxTest::TestSuite {
  /// @brief Number of samples
  int num_samples_;

  /// @brief Number of items in the population
  int num_items_;

  /// @brief Sampler under test
  WeightedReservoirSampler<int>* sampler_;

 public:
  void setUp() {
    num_samples_ = 100;
    num_items_ = 10000;

    // create a synthetic dataset
    sampler_ = new WeightedReservoirSampler<int>(num_samples_);
    for (int i = 1; i <= num_items_; ++i)
      sampler_->consider_sample(i, numeric::random::uniform());
  }

  void tearDown() {
    delete sampler_;
  }

  /// @brief If no items are considered, no samples should be returned
  void test_no_samples_considered() {
    WeightedReservoirSampler<int> sampler(num_samples_);
    vector1<int> selected;
    sampler.samples(&selected);
    TS_ASSERT_EQUALS(0, selected.size());
  }

  /// @brief If no items are considered, no samples should be returned
  void test_more_items_than_samples() {
    vector1<int> selected;
    sampler_->samples(&selected);
    TS_ASSERT_EQUALS(num_samples_, selected.size());
  }
};

}  // anonymous namespace
