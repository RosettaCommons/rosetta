// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/star/Extender.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_STAR_EXTENDER_HH_
#define PROTOCOLS_STAR_EXTENDER_HH_

// External headers
#include <boost/utility.hpp>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

namespace protocols {
namespace star {

class Extender : private boost::noncopyable {
 public:
  Extender(core::sequence::SequenceAlignmentCOP alignment, int num_residues);

  /// @detail Sets unaligned residues' torsion angles to their extended values
  /// and idealizes bond lengths and angles. The placement of unaligned residues
  /// depends on the number of residues separating the adjacent aligned regions.
  ///
  /// Given aligned regions a_1 and a_2 with (a_1.start < a_2.start):
  ///
  /// If (a_2.start - a_1.stop) <= -abinitio:star:short_loop_len, the unaligned
  /// residues are grown off a_1.stop. Otherwise, a stochastically selected
  /// portion of the residues are grown off a_1.stop and the remainder grown off
  /// a_2.start. In both cases, the method is responsible for noting the interior
  /// cutpoints it selected.
  void extend_unaligned(core::pose::Pose* pose);

  /// @brief Returns the unaligned regions in increasing order of start position
  protocols::loops::LoopsCOP unaligned() const {
    return unaligned_;
  }

  /// @brief Returns the aligned regions in increasing order of start position
  protocols::loops::LoopsCOP aligned() const {
    return aligned_;
  }

  /// @brief Returns a set of suggested cutpoints for star fold tree construction
  /// based on information gathered during extend_unaligned().
  const utility::vector1<int>& cutpoints() const {
    return cutpoints_;
  }

  /// @brief Updates predicted secondary structure, improving cutpoint selection
  void set_secondary_structure(core::fragment::SecondaryStructureCOP pred_ss) {
    pred_ss_ = pred_ss;
  }

 protected:
  /// @brief Helper method for settings the torsion angles of a contiguous set
  /// of residues to their extended values and idealizing bond length and angle.
  void extend_section(int start, int stop, core::pose::Pose* pose) const;

  /// @detail Selects a cutpoint on the closed interval [start, stop] using
  /// weighted reservoir sampling. Each residue's weight is proportional to its
  /// likelihood of being a loop.
  int choose_cutpoint(int start, int stop) const;

  /// @detail Replaces pose's fold tree with a simple fold tree of num_residues.
  void simple_fold_tree(core::pose::Pose* pose) const;

 private:
  core::sequence::SequenceAlignmentCOP alignment_;
  protocols::loops::LoopsOP unaligned_;
  protocols::loops::LoopsOP aligned_;
  int num_residues_;

  core::fragment::SecondaryStructureCOP pred_ss_;
  utility::vector1<int> cutpoints_;
};

}  // namespace star
}  // namespace protocols

#endif  // PROTOCOLS_ABINITIO_STAR_EXTENDER_HH_
