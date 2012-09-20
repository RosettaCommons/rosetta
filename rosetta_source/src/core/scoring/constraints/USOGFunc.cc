// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/constraints/USOGFunc.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Package headers
#include <core/scoring/constraints/USOGFunc.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/util.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <algorithm>
#include <cmath>
#include <iostream>

namespace core {
namespace scoring {
namespace constraints {

const core::Real USOGFunc::kDerivativeWindow = 1e-6;
const core::Real USOGFunc::kMinGaussianScore = 1e-8;

USOGFunc::USOGFunc(const utility::vector1<core::Real>& means,
                   const utility::vector1<core::Real>& std_devs,
                   const utility::vector1<core::Real>& weights)
    : means_(means), std_devs_(std_devs), weights_(weights) {
  if (means_.size() != std_devs_.size() || means_.size() != weights_.size()) {
    utility_exit_with_message("Unequal number of means, std_devs, weights");
  }
}

FuncOP USOGFunc::clone() const {
  return new USOGFunc(*this);
}

core::Real USOGFunc::func(const core::Real x) const {
  return -std::log(gaussianScore(x));
}

core::Real USOGFunc::dfunc(const core::Real x) const {
  core::Real w = kDerivativeWindow;
  return (func(x + w) - func(x - w)) / (2 * w);
}

core::Real USOGFunc::gaussianScore(const core::Real x) const {
  core::Real sum = 0;
  for (core::Size i = 1; i <= numGaussians(); ++i) {
    sum += dgaussian(x, means_[i], std_devs_[i], weights_[i]);
  }

  return std::max(sum, kMinGaussianScore);
}

core::Size USOGFunc::numGaussians() const {
  return means_.size();
}

void USOGFunc::show_definition(std::ostream& out) const {
  out << "USOGFUNC " << numGaussians();
  for (core::Size i = 1; i <= numGaussians(); ++i) {
    out << " " << means_[i] << " " << std_devs_[i] << " " << weights_[i];
  }
}

void USOGFunc::resetInstance() {
  means_.clear();
  weights_.clear();
  std_devs_.clear();
}

// Format: USOGFunc <num_gaussians> <mean_1> <std_dev1> <weight_1> ... <mean_n> <std_dev_n> <weight_n>
void USOGFunc::read_data(std::istream& in) {
  resetInstance();

  core::Size num_gaussians = 0;
  in >> num_gaussians;

  for (core::Size i = 1; i <= num_gaussians; ++i) {
    core::Real m = readValueOrDie(in);
    core::Real s = readValueOrDie(in);
    core::Real w = readValueOrDie(in);

    assert(s > 0);
    assert(w > 0);

    means_.push_back(m);
    weights_.push_back(w);
    std_devs_.push_back(s);
  }
}

// Utility functions
core::Real readValueOrDie(std::istream& in) {
  core::Real x;
  in >> x;

  if (in.fail()) {
    utility_exit_with_message("Failed to read floating point value ");
  }

  return x;
}

}  // namespace constraints
}  // namespace scoring
}  // namespace core
