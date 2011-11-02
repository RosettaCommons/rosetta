// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NonlocalAbinitio.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_NONLOCALABINITIO_HH_
#define PROTOCOLS_NONLOCAL_NONLOCALABINITIO_HH_

// Unit header
#include <protocols/nonlocal/NonlocalAbinitio.fwd.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/fragment/FragSet.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/sequence/SequenceAlignment.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/nonlocal/TreeBuilder.fwd.hh>

#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>


namespace protocols {
namespace nonlocal {

class NonlocalAbinitio : public protocols::moves::Mover {
 public:
  NonlocalAbinitio();

  // -- mover -- //
  std::string get_name() const;
  void apply(core::pose::Pose& pose);

  // -- jd2 -- //
  protocols::moves::MoverOP clone() const;
  protocols::moves::MoverOP fresh_instance() const;

private:
  /// @brief Returns a pointer to the large fragment library
  core::fragment::FragSetOP fragments_large() const;

  /// @brief Returns a pointer to the small fragment library
  core::fragment::FragSetOP fragments_small() const;

  /// @brief Identify aligned / unaligned regions by scanning the alignment.
  /// Limit the lengths of these regions to enhance conformational sampling.
  void identify_chunks(const core::sequence::SequenceAlignment& alignment,
                       protocols::loops::Loops* chunks) const;

  /// @brief Estimates missing backbone density
  void build_partial_model(core::pose::Pose* pose) const;

  /// @brief Closes any remaining loops
  void loop_closure(core::pose::Pose* pose) const;

  /// @brief Defines the kinematics of the system
  TreeBuilderOP make_fold_tree(const protocols::loops::Loops& regions, core::pose::Pose* pose) const;

  /// @brief Defines restrictions on degree of freedom modifications
  core::kinematics::MoveMapOP make_movemap(const core::kinematics::FoldTree& tree) const;

  /// @brief Full-atom refinement
  void refine(core::pose::Pose* pose) const;

  // -- members -- //
  core::fragment::FragSetOP fragments_lg_;
  core::fragment::FragSetOP fragments_sm_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_NONLOCALABINITIO_HH_
