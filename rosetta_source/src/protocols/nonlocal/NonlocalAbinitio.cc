// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NonlocalAbinitio.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/NonlocalAbinitio.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/jd2/ThreadingJob.hh>
#include <protocols/loops/LoopRelaxThreadingMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <protocols/relax/FastRelax.hh>

// Package headers
#include <protocols/nonlocal/BrokenFold.hh>
#include <protocols/nonlocal/TreeBuilder.hh>
#include <protocols/nonlocal/TreeBuilderFactory.hh>
#include <protocols/nonlocal/util.hh>

namespace protocols {
namespace nonlocal {

static basic::Tracer TR("protocols.nonlocal.NonlocalAbinitio");

NonlocalAbinitio::NonlocalAbinitio() {
  using core::fragment::FragmentIO;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  FragmentIO io;
  fragments_lg_ = io.read_data(option[in::file::frag9]());
  fragments_sm_ = io.read_data(option[in::file::frag3]());
}

void NonlocalAbinitio::apply(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using protocols::jd2::ThreadingJob;
  using protocols::loops::Loops;
  using protocols::moves::MoverOP;

  ThreadingJob const * const job = current_job();

  // Estimate missing backbone density by performing fragment insertion and
  // loop closure using a simple fold tree
  estimate_missing_density(&pose);

  Loops chunks;
  identify_chunks(job->alignment(), pose.total_residue(), &chunks);
  TreeBuilderOP builder = make_fold_tree(chunks, &pose);

  // Define the kinematics of the system
  emit_intermediate(pose, "nla_pre_abinitio.pdb");
  MoverOP mover = new BrokenFold(fragments_large(), fragments_small(), make_movemap(pose.fold_tree()));
  mover->apply(pose);
  emit_intermediate(pose, "nla_post_abinitio.pdb");

  // Revert any modifications to the pose that TreeBuilder introduced
  builder->tear_down(&pose);
  refine(&pose);
}

TreeBuilderOP NonlocalAbinitio::make_fold_tree(const protocols::loops::Loops& regions,
                                               core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  assert(pose);
  TR << "Regions: " << regions << std::endl;

  TreeBuilderOP builder = TreeBuilderFactory::get_builder(option[OptionKeys::nonlocal::builder]());
  builder->set_up(regions, pose);

  TR << pose->fold_tree() << std::endl;
  return builder;
}

core::kinematics::MoveMapOP NonlocalAbinitio::make_movemap(const core::kinematics::FoldTree& tree) const {
  using core::Size;
  using core::kinematics::MoveMap;
  using core::kinematics::MoveMapOP;

  MoveMapOP movable = new MoveMap();
  movable->set_bb(true);
  movable->set_chi(true);
  movable->set_jump(true);

  // Prevent modification to cutpoint residues and their immediate neighbors
  for (Size i = 1; i <= tree.nres(); ++i) {
    if (tree.is_cutpoint(i)) {
      Size start = std::max(i-2, (Size) 1);
      Size stop  = std::min(i+2, tree.nres());
      for ( Size j = start; j <= stop; ++j ) {
        movable->set_bb(j, false);
      }
    }
  }
  return movable;
}

void NonlocalAbinitio::estimate_missing_density(core::pose::Pose* pose) const {
  assert(pose);

  protocols::loops::LoopRelaxThreadingMover closure;
  closure.setup();
  closure.apply(*pose);

  core::util::switch_to_residue_type_set(*pose, core::chemical::CENTROID);
}

void NonlocalAbinitio::refine(core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Real;
  using core::scoring::ScoreFunctionOP;
  using core::scoring::ScoreFunctionFactory;
  using protocols::relax::FastRelax;

  assert(pose);

  // Relax into constraints according to the search strategy in use
  Real base_wt = option[OptionKeys::constraints::cst_fa_weight]();
  Real constraint_wt = base_wt / 2;

  // Optionally relax
  if (option[OptionKeys::abinitio::relax]()) {
    std::string weights("score12_full");
    if (option[OptionKeys::score::weights].user()) {
      weights = option[score::weights]();
    }

    ScoreFunctionOP score = ScoreFunctionFactory::create_score_function(weights);
    if (option[OptionKeys::score::patch].user()) {
      score->apply_patch_from_file(option[score::patch]());
    }

    // relax with constraints when possible
    if (option[OptionKeys::in::file::rdc].user()) {
      score->set_weight(core::scoring::rdc, option[OptionKeys::nonlocal::rdc_weight]());
    }

    if (option[constraints::cst_file].user()) {
      score->set_weight(core::scoring::atom_pair_constraint, constraint_wt);
      score->set_weight(core::scoring::coordinate_constraint, option[OptionKeys::constraints::cst_weight]());
    }

    FastRelax relax(score);
    relax.apply(*pose);
  }
}

void NonlocalAbinitio::identify_chunks(const core::sequence::SequenceAlignment& alignment,
                                       const core::Size num_residues,
                                       protocols::loops::Loops* chunks) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;
  using protocols::loops::Loops;

  Size min_chunk_sz = fragments_lg_->max_frag_length();
  Size max_chunk_sz = option[OptionKeys::nonlocal::max_chunk_size]();

  Loops aligned, unaligned;
  find_regions_with_minimum_size(alignment, min_chunk_sz, &aligned, &unaligned);
  limit_chunk_size(min_chunk_sz, max_chunk_sz, &aligned);
  limit_chunk_size(min_chunk_sz, max_chunk_sz, &unaligned);

  /// TODO(cmiles) is there anything in alignment that is a proxy for num_residues?
  *chunks = combine_and_trim(min_chunk_sz, num_residues, aligned, unaligned);
  chunks->sequential_order();
}

std::string NonlocalAbinitio::get_name() const {
  return "NonlocalAbinitio";
}

protocols::moves::MoverOP NonlocalAbinitio::clone() const {
  return new NonlocalAbinitio(*this);
}

protocols::moves::MoverOP NonlocalAbinitio::fresh_instance() const {
  return new NonlocalAbinitio();
}

core::fragment::FragSetOP NonlocalAbinitio::fragments_large() const {
  return fragments_lg_;
}

core::fragment::FragSetOP NonlocalAbinitio::fragments_small() const {
  return fragments_sm_;
}

}  // namespace nonlocal
}  // namespace protocols
