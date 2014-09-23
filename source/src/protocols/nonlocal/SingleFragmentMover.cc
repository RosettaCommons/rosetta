// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/SingleFragmentMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/SingleFragmentMover.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>

// Package headers
#include <protocols/nonlocal/Chunk.hh>
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>

typedef protocols::moves::Mover Parent;

namespace protocols {
namespace nonlocal {

static thread_local basic::Tracer TR( "protocols.nonlocal.SingleFragmentMover" );


SingleFragmentMover::SingleFragmentMover() : Parent("SingleFragmentMover") {}

SingleFragmentMover::SingleFragmentMover(const core::fragment::FragSetOP& fragments,
                                         const core::kinematics::MoveMapOP& movable)
    : Parent("SingleFragmentMover") {
  initialize(fragments, movable, PolicyFactory::get_policy("uniform", fragments));
}

SingleFragmentMover::SingleFragmentMover(const core::fragment::FragSetOP& fragments,
                                         const core::kinematics::MoveMapOP& movable,
                                         const PolicyOP& policy)
    : Parent("SingleFragmentMover") {
  initialize(fragments, movable, policy);
}

void SingleFragmentMover::initialize(const core::fragment::FragSetOP& fragments,
                                     const core::kinematics::MoveMapOP& movable,
                                     const PolicyOP& policy) {
  assert(fragments);
  assert(movable);
  assert(policy);

  // Initialize member variables
  fragments_ = fragments;
  movable_ = movable;
  policy_ = policy;

  // Create a position-indexable lookup for core::fragment::Frame's
  initialize_library();
}

void SingleFragmentMover::apply(core::pose::Pose& pose) {
  using core::Size;
  using core::fragment::FragDataCOP;
  using core::fragment::Frame;
  using core::kinematics::FoldTree;

  // ensure that preconditions on the input have been met
  assert(pose.total_residue() > 0);
  assert(pose.fold_tree().check_fold_tree());
  assert(valid());

  // determine whether the pose is in fullatom mode. if so, warn the user and
  // convert it to centroid mode automatically.
  // bool was_fullatom = // Unused variable causes warning.
	to_centroid(&pose);

  // reuse <chunks_> when possible
  const FoldTree& current_tree = pose.fold_tree();
  if (!previous_tree_ || *previous_tree_ != current_tree) {
    chunks_.clear();
    probs_.clear();
    initialize_chunks(current_tree);
    previous_tree_ = FoldTreeOP( new core::kinematics::FoldTree(current_tree) );
  }

  // randomly select the insertion position
  const Chunk* chunk = random_chunk();
  if (!chunk) {
    TR.Warning << "No move possible-- 0 chunks" << std::endl;
    return;
  }

  Size insertion_pos = chunk->choose();

	while (!movable_->get_bb(insertion_pos))
		insertion_pos = chunk->choose();

  // delegate responsibility for choosing the fragment to the policy
  const Frame& frame = library_[insertion_pos];
  Size fragment_pos = policy_->choose(frame, pose);
  frame.apply(*movable_, fragment_pos, pose);
  TR.Debug << "Inserted fragment at position " << insertion_pos << std::endl;
}

std::string SingleFragmentMover::get_name() const {
  return "SingleFragmentMover";
}

protocols::moves::MoverOP SingleFragmentMover::clone() const {
  return protocols::moves::MoverOP( new SingleFragmentMover(*this) );
}

protocols::moves::MoverOP SingleFragmentMover::fresh_instance() const {
  return protocols::moves::MoverOP( new SingleFragmentMover() );
}

void SingleFragmentMover::parse_my_tag(const utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	const protocols::filters::Filters_map&,
	const protocols::moves::Movers_map&,
	const core::pose::Pose& pose) {
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;
  using core::kinematics::MoveMap;
  using core::kinematics::MoveMapOP;
  using std::string;

  // required flags
  // fragments, movemap
  string fragments_file = tag->getOption<string>("fragments");
  FragSetOP fragments = FragmentIO().read_data(fragments_file);

  // initially, all backbone torsions are movable
  MoveMapOP movable( new MoveMap() );
  protocols::rosetta_scripts::parse_movemap(tag, pose, movable, data);

  // optional flags
  // policy -- default => uniform
  string policy_type = tag->getOption<string>("policy", "uniform");
  PolicyOP policy = PolicyFactory::get_policy(policy_type, fragments);

  initialize(fragments, movable, policy);
}

void SingleFragmentMover::parse_def( utility::lua::LuaObject const & def,
	utility::lua::LuaObject const & ,
	utility::lua::LuaObject const & ,
	protocols::moves::MoverCacheSP ) {
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;
  using core::kinematics::MoveMap;
  using core::kinematics::MoveMapOP;
  using std::string;

  // required flags
  // fragments, movemap
  string fragments_file = def["fragments"].to<string>();
  FragSetOP fragments = FragmentIO().read_data(fragments_file);

  // initially, all backbone torsions are movable
  MoveMapOP movable( new MoveMap() );
	movable->set_bb( true );
	movable->set_chi( true );
	movable->set_jump( true );
	if( def["movemap"] )
		protocols::elscripts::parse_movemapdef(def["movemap"], movable);

  // optional flags
  // policy -- default => uniform
  string policy_type = def["policy"] ? def["policy"].to<string>() : "uniform";
  PolicyOP policy = PolicyFactory::get_policy(policy_type, fragments);

  initialize(fragments, movable, policy);
}

bool SingleFragmentMover::valid() const {
  return fragments_ && movable_ && policy_;
}

void SingleFragmentMover::initialize_library() {
  core::fragment::ConstFrameIterator i;
  for (i = fragments_->begin(); i != fragments_->end(); ++i) {
    Size position = (*i)->start();
    library_[position] = **i;
  }
}

void SingleFragmentMover::initialize_chunks(const core::kinematics::FoldTree& tree) {
  using core::Size;
  using core::Real;

  // Assumption: constant-length fragments at each insertion position
  Size fragment_len = fragments_->max_frag_length();

  // Last position to examine
  Size last_pos = fragments_->max_pos() - fragment_len + 1;

  Size p = fragments_->min_pos();
  do {
    Size q = p + 1;
    while (!tree.is_cutpoint(q++)) {}

    // During fold tree construction, it may happen that the distance between
    // the <p> and the next cutpoint is less than fragment length. In this case,
    // the Region::start() exceeds Region::stop(). It's impossible to perform
    // fragment insertions on this region given the current fragment library.
    RegionOP region( new Region(p, q - fragment_len) );
    if (region->increasing()) {
      // Ensure that the chunk is valid before adding it to the list. Mainly, this
      // means that there must be at least 1 movable residue.
      Chunk chunk(region, movable_);

      if (chunk.is_movable()) {
        TR.Debug << "Added chunk: " << region->start() << "-" << region->stop() << std::endl;
        chunks_.push_back(chunk);
      } else {
        TR.Debug << "Skipped chunk: " << region->start() << "-" << region->stop() << ": no movable positions" << std::endl;
      }
    }

    p = q - 1;
  } while (p < last_pos);

  // Assign probabilities of selection to each Chunk according to its length
  Real N = fragments_->max_pos();
  Real sum = 0;
  for (Size i = 1; i <= chunks_.size(); ++i) {
    Real prob = chunks_[i].length() / N;
    probs_.push_back(prob);
    sum += prob;
  }

  // Normalize probabilities
  for (Size i = 1; i <= probs_.size(); ++i) {
    probs_[i] /= sum;
    TR << "P(c_" << i << ")=" << probs_[i] << std::endl;
  }
}

const Chunk* SingleFragmentMover::random_chunk() const {
  using core::Real;
  using core::Size;
  using utility::vector1;

  assert(chunks_.size() > 0);
  vector1<Real> fitnesses(probs_);

  // convert to a CDF
  Size n = fitnesses.size();
  for (Size i = 2; i <= n; ++i)
    fitnesses[i] += fitnesses[i-1];

  // random number from 0 to f[n] (inclusive)
  Real x = numeric::random::uniform();

  // this could be done more efficiently with binary search
  for (Size i = 2; i <= n; ++i) {
    if (fitnesses[i-1] < x && x <= fitnesses[i])
      return &chunks_[i];
  }
  return &chunks_[1];
}

bool SingleFragmentMover::to_centroid(core::pose::Pose* pose) const {
  assert(pose);
  if (pose->is_fullatom()) {
    TR.Warning << "Input pose is full atom (centroid required)" << std::endl;
    TR.Warning << "Performing implicit conversion..." << std::endl;
    core::util::switch_to_residue_type_set(*pose, "centroid");
    return true;
  }
  return false;
}

}  // namespace nonlocal
}  // namespace protocols
