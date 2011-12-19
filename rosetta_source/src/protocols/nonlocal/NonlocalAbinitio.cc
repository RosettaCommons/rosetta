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
// AUTO-REMOVED #include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Project headers
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <core/fragment/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/sequence/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
// AUTO-REMOVED #include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/comparative_modeling/LoopRelaxThreadingMover.hh>
#include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <protocols/loops/util.hh>
#include <protocols/relax/FastRelax.hh>

// Package headers
#include <protocols/nonlocal/BrokenFold.hh>
#include <protocols/nonlocal/TreeBuilder.hh>
#include <protocols/nonlocal/TreeBuilderFactory.hh>
#include <protocols/nonlocal/util.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>



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
  using protocols::comparative_modeling::ThreadingJob;
  using protocols::loops::Loops;

  ThreadingJob const * const job = current_job();

  // Estimate missing backbone density
  build_partial_model(&pose);

  // Decompose the structure into chunks by scanning the alignment.
  // Explicitly limit chunk length to enhance conformational sampling.
  Loops chunks;
  identify_chunks(job->alignment(), &chunks);

  // Define the kinematics
  TreeBuilderOP builder = make_fold_tree(chunks, &pose);

  BrokenFold mover(fragments_large(),
                   fragments_small(),
                   make_movemap(pose.fold_tree()));

  mover.apply(pose);

  // Revert any modifications to the pose that TreeBuilder introduced
  builder->tear_down(&pose);

  // Full-atom refinement
  loop_closure(&pose);
  refine(&pose);
}

TreeBuilderOP NonlocalAbinitio::make_fold_tree(const protocols::loops::Loops& regions,
                                               core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  assert(pose);

  TreeBuilderOP builder = TreeBuilderFactory::get_builder(option[OptionKeys::nonlocal::builder]());
  builder->set_up(regions, pose);

  TR << pose->fold_tree() << std::endl;
  return builder;
}

core::kinematics::MoveMapOP NonlocalAbinitio::make_movemap(const core::kinematics::FoldTree& tree) const {
	core::kinematics::MoveMapOP movable = new core::kinematics::MoveMap();
  movable->set_bb(true);
  movable->set_jump(true);
  return movable;
}

void NonlocalAbinitio::build_partial_model(core::pose::Pose* pose) const {
  assert(pose);

  protocols::comparative_modeling::LoopRelaxThreadingMover closure;
  closure.setup();
  closure.apply(*pose);

  core::util::switch_to_residue_type_set(*pose, core::chemical::CENTROID);
}

void NonlocalAbinitio::loop_closure(core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::kinematics::FoldTree;
  using protocols::comparative_modeling::LoopRelaxMover;
  using protocols::loops::LoopsOP;
  assert(pose);

  // An empty Loops selection notifies LoopRelaxMover that it is responsible
  // for choosing the breaks to close automatically.
  LoopsOP empty = new protocols::loops::Loops();
  LoopRelaxMover closure;
  closure.remodel("quick_ccd");
  closure.intermedrelax("no");
  closure.refine("no");
  closure.relax("no");
  closure.loops(empty);

  // 3-mers
  utility::vector1<core::fragment::FragSetOP> fragments;
  fragments.push_back(fragments_small());
  closure.frag_libs(fragments);

  // Use atom pair constraints when available
  closure.cmd_line_csts(option[constraints::cst_fa_file].user());

  FoldTree tree(pose->total_residue());
  pose->fold_tree(tree);
  closure.apply(*pose);
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
  const Real base_wt = option[OptionKeys::constraints::cst_fa_weight]();
  const Real constraint_wt = base_wt / 2;

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
                                       protocols::loops::Loops* chunks) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;
  using protocols::loops::Loops;

  const Size min_chunk_sz = fragments_lg_->max_frag_length();
  const Size max_chunk_sz = option[OptionKeys::nonlocal::max_chunk_size]();

  protocols::loops::LoopsOP aligned = new protocols::loops::Loops();
  protocols::loops::LoopsOP unaligned = new protocols::loops::Loops();
  
  find_regions_with_minimum_size(alignment, min_chunk_sz, aligned, unaligned);
  limit_chunk_size(min_chunk_sz, max_chunk_sz, aligned);
  limit_chunk_size(min_chunk_sz, max_chunk_sz, unaligned);

  *chunks = combine_and_trim(min_chunk_sz, alignment.length(), aligned, unaligned);
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
