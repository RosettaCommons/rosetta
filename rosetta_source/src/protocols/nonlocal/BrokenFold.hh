// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BrokenFold.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_BROKENFOLD_HH_
#define PROTOCOLS_NONLOCAL_BROKENFOLD_HH_

// Unit header
#include <protocols/nonlocal/BrokenFold.fwd.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace nonlocal {

/// @class Broken-chain protocol for NonlocalAbinitio. Movers are
/// enqueued in the constructor and executed in order in the base
/// class. If specialization is required, simply override apply().
class BrokenFold : public protocols::moves::Mover {
 public:
  BrokenFold(core::fragment::FragSetOP fragments_lg,
             core::fragment::FragSetOP fragments_sm,
             core::kinematics::MoveMapOP movable);

  /// @brief Performs alternating rounds of jump minimization and
  /// broken-chain folding
  void apply(core::pose::Pose& pose);

  /// @brief Returns the name of this mover
  std::string get_name() const;

 private:
  /// @brief Returns the number of cycles for stage i
  core::Size cycles(int stage);

  /// @brief Configures the score function to be minimized for the ith stage
  core::scoring::ScoreFunctionOP score_function(int stage, int num_residues) const;

  /// @brief Writes the stage header to <out>
  void show_stage_header(core::Size stage_num, std::ostream& out) const;

  /// @brief Collection of movers to be applied in order
  utility::vector1<protocols::moves::MoverOP> movers_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_BROKENFOLD_HH_
