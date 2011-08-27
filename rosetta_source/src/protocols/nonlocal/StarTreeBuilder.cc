// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/StarTreeBuilder.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/StarTreeBuilder.hh>

// C/C++ headers
#include <iostream>
#include <utility>

// External headers
#include <boost/format.hpp>

// Objexx headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/WeightedReservoirSampler.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/StructuralConservationStore.hh>

// Package headers
#include <protocols/nonlocal/CutFinder.hh>
#include <protocols/nonlocal/NLGrouping.hh>

namespace protocols {
namespace nonlocal {

static basic::Tracer TR("protocols.nonlocal.StarTreeBuilder");

StarTreeBuilder::StarTreeBuilder() {
  virtual_residue(-1);
}

// TODO(cmiles) Method getting a little long, consider refactoring
//
// Assumes <grouping> is sorted
void StarTreeBuilder::set_up(const NLGrouping& grouping, core::pose::Pose* pose) {
  using core::Size;
  using core::Real;
  using core::fragment::SecondaryStructure;
  using core::fragment::SecondaryStructureOP;
  using numeric::random::WeightedReservoirSampler;
  using utility::vector1;
  assert(pose->fold_tree().check_fold_tree());

  // Number of residues before addition of virtual residue
  Size num_residues = pose->total_residue();

  // Add a virtual residue at the center of mass, updating the FoldTree
  numeric::xyzVector<Real> center;
  grouping.center_of_mass(&center);
  core::pose::addVirtualResAsRoot(center, *pose);
  TR.Debug << "center of mass(aligned region): "
           << boost::format("(%1%, %2%, %3%)") % center.x() % center.y() % center.z()
           << std::endl;

  // Initialize member variable <virtual_res_> with the index of the newly added
  // virtual residue. Subsequently, <virtual_res_> can serve as a proxy for
  // <num_residues>
  virtual_residue(pose->total_residue());
  core::kinematics::FoldTree tree(pose->fold_tree());
  TR.Debug << "starting fold tree: " << tree << std::endl;

  bool has_conservation = pose->data().has(core::pose::datacache::CacheableDataType::STRUCTURAL_CONSERVATION);
  TR << "Structural conservation available: " << has_conservation << std::endl;

  utility::vector1<std::pair<int, int> > jumps;
  for (Size i = 1; i <= grouping.num_groups(); ++i) {
    const NLFragmentGroup& fragment = grouping.groups(i);

    // When available, use structural conservation to inform jump selection.
    // Otherwise, choose the midpoint of the NLFragmentGroup.
    Size anchor_position = has_conservation ?
        choose_conserved_position(fragment, *pose) :
        (fragment.start() + fragment.stop()) / 2;

    // virtual residue => anchor position
    jumps.push_back(std::make_pair(virtual_residue(), anchor_position));
    TR << "Added jump: " << virtual_residue() << " => " << anchor_position << std::endl;
  }

  // TODO(cmiles) Improve inter-jump cutpoint selection
  //
  // Create a cut between pairs of jumps. After the loop, add a cut between the
  // final jump and the end of the chain.
  SecondaryStructureOP ss = new SecondaryStructure(*pose);
  utility::vector1<int> cuts;
  for (Size i = 2; i <= jumps.size(); ++i) {
    Size cutpoint = CutFinder::choose_cutpoint(jumps[i-1].second + 1, jumps[i].second - 1, ss);
    cuts.push_back(cutpoint);
    TR.Debug << "added cut at " << cutpoint << std::endl;
  }

  // Remember to include the original cutpoint at the end of the chain
  // (before the virtual residue)
  cuts.push_back(num_residues);

  // Construct the star fold tree from the set of jumps and cuts above.
  // Reorder the resulting fold tree so that <virtual_res> is the root.
  ObjexxFCL::FArray2D_int ft_jumps(2, jumps.size());
  for (Size i = 1; i <= jumps.size(); ++i) {
    ft_jumps(1, i) = std::min(jumps[i].first, jumps[i].second);
    ft_jumps(2, i) = std::max(jumps[i].first, jumps[i].second);
  }

  ObjexxFCL::FArray1D_int ft_cuts(cuts.size());
  for (Size i = 1; i <= cuts.size(); ++i) {
    ft_cuts(i) = cuts[i];
  }

  // Try to build the star fold tree from jumps and cuts
  bool status = tree.tree_from_jumps_and_cuts(virtual_residue(),   // nres_in
                                              jumps.size(),        // num_jump_in
                                              ft_jumps,            // jump_point
                                              ft_cuts,             // cuts
                                              virtual_residue());  // root
  if (!status)
    utility_exit_with_message("StarTreeBuilder: failed to build fold tree from cuts and jumps");

  // Update the pose's fold tree
  pose->fold_tree(tree);
}

core::Size StarTreeBuilder::choose_conserved_position(const NLFragmentGroup& fragment,
                                                      const core::pose::Pose& pose) {
  utility::vector1<core::Size> positions;
  utility::vector1<core::Real> weights;
  for (core::Size j = fragment.start(); j <= fragment.stop(); ++j) {
    positions.push_back(j);
    weights.push_back(core::pose::structural_conservation(pose, j));
  }

  // Perform weighted selection to choose (a single) jump position
  utility::vector1<core::Size> samples;
  numeric::random::WeightedReservoirSampler<Size> sampler(1);
  for (core::Size j = 1; j <= positions.size(); ++j) {
    sampler.consider_sample(positions[j], weights[j]);
  }
  sampler.samples(&samples);
  return samples[1];
}

void StarTreeBuilder::tear_down(core::pose::Pose* pose) {
  if (virtual_res_ == -1) {
    TR.Warning << "Attempt to tear_down() unitialized virtual residue" << std::endl;
    return;
  } else {
    pose->conformation().delete_residue_slow(virtual_res_);
    TR << "Removed virtual residue at position " << virtual_res_ << std::endl;
    virtual_residue(-1);
  }
}

int StarTreeBuilder::virtual_residue() const {
  return virtual_res_;
}

void StarTreeBuilder::virtual_residue(int virtual_res) {
  virtual_res_ = virtual_res;
}


}  // namespace nonlocal
}  // namespace protocols
