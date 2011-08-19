// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/distribution.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>

// External headers
#include <boost/math/distributions.hpp>

// Utility headers
#include <devel/init.hh>
#include <numeric/random/DistributionSampler.hh>

using namespace std;
using namespace boost::math;
using numeric::random::DistributionSampler;

#define NUM_SAMPLES 100

void sample_normal() {
  normal dist;
  DistributionSampler<normal> sampler(dist);

  cout << "Normal:" << endl;
  for (int i = 1; i <= NUM_SAMPLES; ++i) {
    cout << sampler.sample() << endl;
  }
  cout << endl;
}

void sample_exponential() {
  exponential dist;
  DistributionSampler<exponential> sampler(dist);

  cout << "Exponential:" << endl;
  for (int i = 1; i <= NUM_SAMPLES; ++i) {
    cout << sampler.sample() << endl;
  }
  cout << endl;
}

int main(int argc, char* argv[]) {
  devel::init(argc, argv);

  // a few illustrative cases...
  sample_normal();
  sample_exponential();
}
