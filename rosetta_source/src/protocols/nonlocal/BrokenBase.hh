// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BrokenBase.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_BROKENBASE_HH_
#define PROTOCOLS_NONLOCAL_BROKENBASE_HH_

// Unit header
#include <protocols/nonlocal/BrokenBase.fwd.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace nonlocal {

/// @class Abstract base class for broken-chain folding specializations
class BrokenBase : public protocols::moves::Mover {
  typedef protocols::moves::MoverOP MoverOP;

 public:
  /// @brief Modifies the pose
  virtual void apply(core::pose::Pose& pose);

  /// @brief Returns the name of this mover
  virtual std::string get_name() const = 0;

 protected:
  /// @brief Adds the specified mover to the collection
  void add_mover(MoverOP mover);

  /// @brief Adds cutpoint variants to <pose>
  void add_cutpoint_variants(core::pose::Pose* pose) const;

  /// @brief Removes cutpoint variants from <pose>
  void remove_cutpoint_variants(core::pose::Pose* pose) const;

  /// @brief Write intermediate pose to disk
  void emit_intermediate(const core::pose::Pose& pose, core::Size stage_num) const;

  /// @brief Writes the stage header to <out>
  void show_stage_header(core::Size stage_num, std::ostream& out) const;

  /// @brief Allocates, configures, and returns a score function for stage i
  core::scoring::ScoreFunctionOP score_function(int stage, int num_residues) const;

  /// @brief Configures a minimizer based on the current score function
  protocols::moves::MoverOP make_minimizer(core::scoring::ScoreFunctionOP score);

 private:
  /// @brief Collection of movers to be applied in order
  utility::vector1<MoverOP> stage_movers_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_BROKENBASE_HH_
