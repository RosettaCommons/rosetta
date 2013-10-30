// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/random/DistributionSampler.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_NUMERIC_RANDOM_DISTRIBUTIONSAMPLER_HH
#define INCLUDED_NUMERIC_RANDOM_DISTRIBUTIONSAMPLER_HH

// External headers
#include <boost/utility.hpp>

#include <boost/math/distributions.hpp>

// Utility headers
#include <numeric/random/random.hh>

namespace numeric {
namespace random {

template <typename T>
class DistributionSampler : boost::noncopyable {
 public:
  /// @brief Creates a new instance that allows samples to be drawn randomly
  /// from <distribution>
  explicit DistributionSampler(T& distribution) : distribution_(distribution) {}

  /// @brief Returns a random value drawn from the distribution
  /// @detail A general method to generate random numbers from an arbitrary
  /// distribution that has a cdf without jumps is to use the inverse function
  /// to the cdf: G(y)=F^{-1}(y). If u(1), ..., u(n) are random numbers from the
  /// uniform on (0,1) distribution then G(u(1)), ..., G(u(n)) is a random
  /// sample from the distribution with cdf F(x).
  double sample() {
    return boost::math::quantile(distribution(), numeric::random::uniform());
  }

  /// @brief Returns the distribution from which this sampler generates values
  const T& distribution() const {
    return distribution_;
  }

 private:
  /// @brief Distribution from which to draw samples
  T distribution_;
};

}  // namespace random
}  // namespace numeric

#endif  // INCLUDED_NUMERIC_RANDOM_DISTRIBUTIONSAMPLER_HH
