// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/SmoothPolicy.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_SMOOTHPOLICY_HH_
#define PROTOCOLS_NONLOCAL_SMOOTHPOLICY_HH_

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/basic_moves/GunnCost.hh>

// Package headers
#include <protocols/nonlocal/Policy.hh>

namespace protocols {
namespace nonlocal {

/// @class Implements the Policy interface. Given a Frame, chooses the fragment
/// that, when applied to the pose, minimizes total distortion ("smooth move").
/// Precomputes fragment endpoint positions in the constructor to improve
/// performance.
class SmoothPolicy : public Policy {
  typedef core::fragment::FragSetCOP FragSetCOP;
  typedef protocols::basic_moves::GunnCost GunnCost;

 public:

  /// @class Simple container that associates fragment indices with Gunn scores
  class Candidate {
    typedef core::Real Real;
    typedef core::Size Size;

   public:
    Candidate(Real score, Size fragment_num)
        : score_(score), fragment_num_(fragment_num) {}

    /// @brief Returns the candidate's score
    Real score() const {
      return score_;
    }

    /// @brief Returns the candidate's position within the Frame
    Size fragment_num() const {
      return fragment_num_;
    }

   private:
    Real score_;
    Size fragment_num_;
  };


  /// @brief Precomputes the endpoints of each fragment, storing the result in a
  /// two-level map indexed by frame- and fragment-number, respectively.
  explicit SmoothPolicy(FragSetCOP fragments);

  /// @brief Given the current state of <pose>, selects the fragment in <frame>
  /// that minimizes overall distortion
  virtual core::Size choose(const core::fragment::Frame& frame,
                            const core::pose::Pose&);


 private:
  GunnCost scorer_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_SMOOTHPOLICY_HH_
