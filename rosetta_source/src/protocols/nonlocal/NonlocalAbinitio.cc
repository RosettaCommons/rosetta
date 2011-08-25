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
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/fragment/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/loops/LoopRelaxMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <protocols/relax/FastRelax.hh>

// Package headers
#include <protocols/nonlocal/BoundaryFinder.hh>
#include <protocols/nonlocal/BrokenBase.hh>
#include <protocols/nonlocal/BrokenRefine.hh>
#include <protocols/nonlocal/NLFragment.hh>
#include <protocols/nonlocal/NLFragmentGroup.hh>
#include <protocols/nonlocal/NLGrouping.hh>
#include <protocols/nonlocal/TreeBuilder.hh>
#include <protocols/nonlocal/TreeBuilderFactory.hh>
#include <protocols/nonlocal/NonlocalAbinitioReader.hh>
#include <protocols/nonlocal/util.hh>

namespace protocols {
namespace nonlocal {

// -- Convenience types -- //
typedef utility::vector1<NLGrouping> NonlocalGroupings;

static basic::Tracer TR("protocols.nonlocal.NonlocalAbinitio");
static numeric::random::RandomGenerator RG(764443637);

void NonlocalAbinitio::check_required_options() const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  std::string prefix = "Failed to specify required option ";

  if (!option[in::file::frag3].user())
    utility_exit_with_message(prefix + "in:file:frag3");
  if (!option[in::file::frag9].user())
    utility_exit_with_message(prefix + "in:file:frag9");
}

NonlocalAbinitio::NonlocalAbinitio() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // Ensure that the required options have been specified
  if (!option[OptionKeys::nonlocal::moves].user())
    utility_exit_with_message("Failed to specify required option nonlocal:moves");

  // Deserialize the pairings file
  NonlocalGroupings groupings;
  NonlocalAbinitioReader::read(option[OptionKeys::nonlocal::moves](), &groupings);

  initialize(groupings);
}

NonlocalAbinitio::NonlocalAbinitio(const NonlocalGroupings& groupings) {
  initialize(groupings);
}

void NonlocalAbinitio::initialize(const NonlocalGroupings& groupings) {
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;
  using core::fragment::SecondaryStructure;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // Ensure that required options have been specified
  check_required_options();

  // Load fragment libraries from file
  FragmentIO io;
  fragments_lg_ = io.read_data(option[in::file::frag9]());
  fragments_sm_ = io.read_data(option[in::file::frag3]());

  // Predict secondary structure from the fragment files
  secondary_struct_ = new SecondaryStructure(*fragments_sm_);

  // Initialize members
  groupings_ = groupings;
}

void NonlocalAbinitio::apply(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;
  using core::kinematics::MoveMap;
  using core::kinematics::MoveMapOP;
  using protocols::moves::MoverOP;

  // Always operate in centroid mode
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);

  // Are we processing a comparative modeling target?
  bool comparative_modeling_input =
      option[in::file::alignment].user() &&
      option[in::file::template_pdb].user() &&
      option[cm::aln_format].user();

  // Randomly select one of the NLGrouping's as the basis for this trajectory
  NLGrouping grouping = groupings_[RG.random_range(1, groupings_.size())];

  TreeBuilderOP builder;
  if (comparative_modeling_input) {
    initial_closure(&pose);
    builder = make_fold_tree(grouping, &pose);
  } else {
    // Insert the torsions from the randomly selected NLGrouping
    for (Size i = 1; i <= grouping.num_groups(); ++i) {
      const NLFragmentGroup& group = grouping.groups(i);

      for (Size j = 1; j <= group.num_entries(); ++j) {
        const NLFragment& entry = group.entries(j);
        Size position = entry.position();
        pose.set_phi(position, entry.phi());
        pose.set_psi(position, entry.psi());
        pose.set_omega(position, entry.omega());
      }
    }

    // Construct the fold tree and orient the non-local fragments
    builder = make_fold_tree(grouping, &pose);
    superimpose(grouping, &pose);
  }

  // If we're operating in RIGID mode, enforce user-defined restrictions on
  // backbone torsion modification
  MoveMapOP movable = new MoveMap();
  prepare_movemap(grouping, pose, movable);

  // Perform fragment-based assembly
  emit_intermediate(pose, "nla_pre_abinitio.pdb");
  MoverOP mover = new BrokenRefine(fragments_small(), movable);
  mover->apply(pose);
  emit_intermediate(pose, "nla_post_abinitio.pdb");

  // At the conclusion of folding, ask the TreeBuilder to revert any modifications
  // to the pose that it introduced (e.g. virtual residues). Afterwards, proceed
  // to close any remaining chainbreaks and (optionally) relax.
  builder->tear_down(&pose);
  final_closure(&pose);
  relax(&pose);

  // Track the provenance of the non-local pairing by inserting an entry into
  // the pose's DataCache. This allows us to determine which starting points
  // (i.e. set of non-local contacts) contribute to low-scoring structures
  core::pose::add_comment(pose, "Nonlocal Provenance", grouping.provenance());
}

TreeBuilderOP NonlocalAbinitio::make_fold_tree(const NLGrouping& grouping,
                                               core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  TreeBuilderOP builder = TreeBuilderFactory::get_builder(option[OptionKeys::nonlocal::builder]());
  builder->build(grouping, pose, secondary_struct_);

  // Ensure that the FoldTree is left in a consistent state
  const core::kinematics::FoldTree& tree = pose->fold_tree();
  assert(tree.check_fold_tree());
  TR << tree << std::endl;

  return builder;
}

void NonlocalAbinitio::prepare_movemap(const NLGrouping& grouping,
                                       const core::pose::Pose& pose,
                                       core::kinematics::MoveMapOP movable) {
  using core::Size;

  movable->set_bb(true);
  movable->set_chi(true);
  movable->set_jump(true);

  // Prevent modification to cutpoint residues and their neighbors
  const core::kinematics::FoldTree& tree = pose.fold_tree();
  for (Size i = 1; i <= pose.total_residue(); ++i) {
    if (tree.is_cutpoint(i)) {
      Size start = std::max(i-2, (Size) 1);
      Size stop  = std::min(i+2, pose.total_residue());
      for ( Size j = start; j <= stop; ++j ) {
        movable->set_bb(j, false);
      }
    }
  }

  // Explicitly flush the tracer's buffer
  movable->show(TR.Debug, pose.total_residue());
  TR.Debug << std::endl;
}

void NonlocalAbinitio::superimpose(const NLGrouping& grouping, core::pose::Pose* pose) const {
  using core::Real;
  using core::Size;
  using numeric::xyzVector;

  for (Size i = 1; i <= grouping.num_groups(); ++i) {
    const NLFragmentGroup& group = grouping.groups(i);

    // Construct stubs from the 3 central CA atoms
    Size central_residue = static_cast<Size>(ceil(group.num_entries() / 2.0));

    const NLFragment& f1 = group.entries(central_residue - 1);
    const NLFragment& f2 = group.entries(central_residue);
    const NLFragment& f3 = group.entries(central_residue + 1);
    xyzVector<Real> fa1(f1.x(), f1.y(), f1.z());
    xyzVector<Real> fa2(f2.x(), f2.y(), f2.z());
    xyzVector<Real> fa3(f3.x(), f3.y(), f3.z());

    Size midpoint_pose = f2.position();
    xyzVector<Real> m1 = pose->residue(f1.position()).xyz("CA");
    xyzVector<Real> m2 = pose->residue(f2.position()).xyz("CA");
    xyzVector<Real> m3 = pose->residue(f3.position()).xyz("CA");

    // Compute the transform
    core::kinematics::Stub x = core::fragment::getxform(m1, m2, m3, fa1, fa2, fa3);

    // Starting at the midpoint of the insertion point in the pose, propagate
    // the change to the left and right until we reach either the end of the
    // chain or a cutpoint
    Size region_start, region_stop;
    BoundaryFinder::boundaries(pose->fold_tree(), midpoint_pose, &region_start, &region_stop);
    core::fragment::xform_pose(*pose, x, region_start, region_stop);
  }
  TR << "Pose reorientation complete!" << std::endl;
}

void NonlocalAbinitio::initial_closure(core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::kinematics::FoldTree;
  using protocols::loops::LoopRelaxMover;
  using protocols::loops::Loops;
  assert(pose);

  // An empty Loops selection notifies LoopRelaxMover that it is responsible
  // for choosing the breaks to close automatically.
  Loops empty;
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

  // LoopRelax constructs and manages its own FoldTree, which may be in
  // conflict with the existing FoldTree. Set the Pose's FoldTree to the
  // desired end state (a simple fold tree).
  FoldTree orig_tree = pose->fold_tree();
  pose->fold_tree(FoldTree(pose->total_residue()));

  // Begin loop closure
  closure.apply(*pose);

  // Restore <pose>'s original fold tree
  pose->fold_tree(orig_tree);
}

void NonlocalAbinitio::final_closure(core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::kinematics::FoldTree;
  using protocols::loops::LoopRelaxMover;
  using protocols::loops::Loops;

  assert(pose);
  pose->fold_tree(FoldTree(pose->total_residue()));

  // Close any remaining chainbreaks
  if (has_chainbreaks(*pose)) {
    // An empty Loops selection notifies LoopRelaxMover that it is responsible
    // for choosing the breaks to close automatically.
    Loops empty;
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

    // LoopRelax constructs and manages its own FoldTree, which may be in
    // conflict with the existing FoldTree. Set the Pose's FoldTree to the
    // desired end state (a simple fold tree).
    pose->fold_tree(FoldTree(pose->total_residue()));

    // Begin loop closure
    closure.apply(*pose);
    emit_intermediate(*pose, "nla_final_closure.pdb");
  }
}

void NonlocalAbinitio::relax(core::pose::Pose* pose) const {
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

std::string NonlocalAbinitio::get_name() const {
  return "NonlocalAbinitio";
}

protocols::moves::MoverOP NonlocalAbinitio::clone() const {
  return new NonlocalAbinitio(*this);
}

protocols::moves::MoverOP NonlocalAbinitio::fresh_instance() const {
  return new NonlocalAbinitio();
}

// -- Accessors -- //

const NonlocalGroupings& NonlocalAbinitio::groupings() const {
  return groupings_;
}

core::fragment::FragSetOP NonlocalAbinitio::fragments_large() const {
  return fragments_lg_;
}

core::fragment::FragSetOP NonlocalAbinitio::fragments_small() const {
  return fragments_sm_;
}

}  // namespace nonlocal
}  // namespace protocols
