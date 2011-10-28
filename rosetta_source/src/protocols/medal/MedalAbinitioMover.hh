// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalAbinitioMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_MEDAL_MEDAL_ABINITIO_MOVER_HH_
#define PROTOCOLS_MEDAL_MEDAL_ABINITIO_MOVER_HH_

// Unit header
#include <protocols/medal/MedalAbinitioMover.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/fragment/FragSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/nonlocal/Policy.fwd.hh>
#include <protocols/nonlocal/StarTreeBuilder.fwd.hh>

// Package headers
#include <protocols/medal/MedalMover.hh>

namespace protocols {
namespace medal {

/// @class Ab initio protocol with medal-like moves and kinematics
class MedalAbinitioMover : public MedalMover {
 public:
  MedalAbinitioMover();
  void apply(core::pose::Pose& pose);

  // -- jd2 -- //
  std::string get_name() const;
  protocols::moves::MoverOP clone() const;
  protocols::moves::MoverOP fresh_instance() const;

 private:
  /// @brief Derive reasonable estimates of residue positions through
  /// fragment insertion with minimal score function
  void estimate_residue_positions(core::pose::Pose* pose) const;

  /// @brief Updates the position of the virtual residue by placing it
  /// at the recomputed center of mass.
  void update_position_vres(const protocols::loops::Loops& chunks,
                            protocols::nonlocal::StarTreeBuilder* builder,
                            core::pose::Pose* pose) const;

	protocols::moves::MoverOP fragment_mover(const core::pose::Pose& pose,
																					 const core::fragment::FragSetOP fragments,
																					 const protocols::nonlocal::PolicyOP policy) const;
};

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_MEDAL_ABINITIO_MOVER_HH_
