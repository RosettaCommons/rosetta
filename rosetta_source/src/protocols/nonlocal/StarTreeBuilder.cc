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
#include <boost/math/distributions/normal.hpp>
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <numeric/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/DistributionSampler.hh>
#include <numeric/random/random.hh>
#include <numeric/random/WeightedReservoirSampler.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/StructuralConservationStore.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/nonlocal/CutFinder.hh>

namespace protocols {
namespace nonlocal {

/// @brief Simple utility method to aid debugging failed fold tree construction
void show_cuts_and_jumps(const utility::vector1<int>& cuts,
                         const utility::vector1<std::pair<int, int> >& jumps) {
  using std::cerr;
  using std::endl;
  using std::pair;

  cerr << "Cutpoints: " << endl;
  for (utility::vector1<int>::const_iterator i = cuts.begin(); i != cuts.end(); ++i) {
    cerr << "\t" << *i << endl;
  }

  cerr << "Jumps: " << endl;
  for (utility::vector1<pair<int, int> >::const_iterator i = jumps.begin(); i != jumps.end(); ++i) {
    const pair<int, int>& jump = *i;
    cerr << "\t" << jump.first << " => " << jump.second << endl;
  }
}

StarTreeBuilder::StarTreeBuilder() : virtual_res_(-1) {}

/// Note: assumes <chunks> are sorted in increasing order of start position
void StarTreeBuilder::set_up(const protocols::loops::Loops& chunks, core::pose::Pose* pose) {
  using core::Size;
  using core::Real;
  using core::fragment::SecondaryStructure;
  using core::fragment::SecondaryStructureOP;
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  using utility::vector1;

  assert(pose);
  assert(chunks.num_loop());

  // Number of residues before addition of virtual residue
  Size num_residues = pose->total_residue();

  // Add a virtual residue at the center of mass (updates the fold tree)
  numeric::xyzVector<Real> center;
  chunks.center_of_mass(*pose, &center);
  core::pose::addVirtualResAsRoot(center, *pose);

  // Initialize member variable <virtual_res_> with the index of the newly added
  // virtual residue. Subsequently, <virtual_res_> can serve as a proxy for
  // <num_residues>
  virtual_res_ = pose->total_residue();
  bool has_conservation = pose->data().has(core::pose::datacache::CacheableDataType::STRUCTURAL_CONSERVATION);

  vector1<std::pair<int, int> > jumps;
  for (Loops::const_iterator i = chunks.begin(); i != chunks.end(); ++i) {
    const Loop& chunk = *i;

    // When available, use structural conservation to inform jump selection.
    // Otherwise, choose the midpoint of the region.
    Size anchor_position = has_conservation ?
        choose_conserved_position(chunk, *pose) :
        choose_unconserved_position(chunk.start(), chunk.stop());

    // virtual residue => anchor position
    jumps.push_back(std::make_pair(virtual_res_, anchor_position));
  }

  // Insert cutpoints between adjacent jumps
  SecondaryStructureOP ss = new SecondaryStructure(*pose);
  vector1<int> cuts;
  for (Size i = 2; i <= jumps.size(); ++i) {
    Size cutpoint = CutFinder::choose_cutpoint(jumps[i-1].second + 1, jumps[i].second - 1, ss);
    cuts.push_back(cutpoint);
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
  core::kinematics::FoldTree tree(pose->fold_tree());
  bool status = tree.tree_from_jumps_and_cuts(virtual_res_,   // nres_in
                                              jumps.size(),   // num_jump_in
                                              ft_jumps,       // jump_point
                                              ft_cuts,        // cuts
                                              virtual_res_);  // root
  if (!status) {
    show_cuts_and_jumps(cuts, jumps);
    utility_exit_with_message("StarTreeBuilder: failed to build fold tree from cuts and jumps");
  }

  // Update the pose's fold tree
  pose->fold_tree(tree);

  // When native is available, compute RMSD over jump residues. Results stored
  // as comments in the pose with format: rmsd_jump_residue_X rmsd
  do_compute_jump_rmsd(pose);
}

void StarTreeBuilder::do_compute_jump_rmsd(core::pose::Pose* model) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Real;
  using core::Size;
  using core::pose::Pose;
  assert(model);

  if (!option[OptionKeys::in::file::native].user())
    return;

  // Compute RMSD of jump residues (jump point +/- 1 residue)
  boost::unordered_map<Size, Real> rmsds;
  Pose native = *core::import_pose::pose_from_pdb(option[OptionKeys::in::file::native]());
  core::scoring::compute_jump_rmsd(native, *model, &rmsds);

  // Write results to <model> as comments
  for (boost::unordered_map<Size, Real>::const_iterator i = rmsds.begin(); i != rmsds.end(); ++i) {
    const Size residue = i->first;
    const Real rmsd = i->second;
    core::pose::add_comment(*model,
                            (boost::format("rmsd_jump_residue_%1%") % residue).str(),
                            (boost::format("%1%") % rmsd).str());
  }
}

// TODO(cmiles) additional testing needed
core::Size StarTreeBuilder::choose_conserved_position(const protocols::loops::Loop& chunk,
                                                      const core::pose::Pose& pose) const {
  using core::Size;
  using utility::vector1;

  vector1<Size> positions;
  vector1<core::Real> weights;
  for (Size j = chunk.start(); j <= chunk.stop(); ++j) {
    positions.push_back(j);
    weights.push_back(core::pose::structural_conservation(pose, j));
  }

  // Perform weighted selection to choose (a single) jump position
  vector1<Size> samples;
  numeric::random::WeightedReservoirSampler<Size> sampler(1);
  for (Size j = 1; j <= positions.size(); ++j) {
    sampler.consider_sample(positions[j], weights[j]);
  }
  sampler.samples(&samples);

  Size position = samples[1];
  return numeric::clamp<Size>(position, chunk.start(), chunk.stop());
}

core::Size StarTreeBuilder::choose_unconserved_position(core::Size start, core::Size stop) const {
  using boost::math::normal;
  using core::Size;
  using numeric::random::DistributionSampler;

  double mu = (stop - start + 1) / 2.0;
  double sigma = 1;
  normal distribution(mu, sigma);
  DistributionSampler<normal> sampler(distribution);

  // Clamp insertion position to closed interval [start, stop]
  Size position = static_cast<Size>(sampler.sample());
  return numeric::clamp<Size>(position, start, stop);;
}

void StarTreeBuilder::tear_down(core::pose::Pose* pose) {
  assert(pose);
  if (virtual_res_ == -1) {
    return;
  } else {
    pose->conformation().delete_residue_slow(virtual_res_);
    virtual_res_ = -1;
  }
}

}  // namespace nonlocal
}  // namespace protocols
