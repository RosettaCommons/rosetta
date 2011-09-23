// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/util.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/util.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/SequenceMapping.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/ThreadingJob.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

namespace protocols {
namespace nonlocal {

static numeric::random::RandomGenerator RG(156144120);
static basic::Tracer TR("protocols.nonlocal.util");

protocols::loops::Loops combine_and_trim(core::Size min_chunk_sz,
                                         core::Size num_residues,
                                         const protocols::loops::Loops& aligned_regions,
                                         const protocols::loops::Loops& unaligned_regions) {
  using protocols::loops::Loops;

  Loops combined;
  for (Loops::const_iterator i = aligned_regions.begin(); i != aligned_regions.end(); ++i)
    combined.push_back(*i);
  for (Loops::const_iterator i = unaligned_regions.begin(); i != unaligned_regions.end(); ++i)
    combined.push_back(*i);

  // trim back the final loop to nres() - min_chunk_sz
  combined[combined.size()].set_stop(num_residues - min_chunk_sz + 1);
  return combined;
}

void limit_chunk_size(core::Size min_chunk_sz,
                      core::Size max_chunk_sz,
                      protocols::loops::Loops* regions) {
  using core::Size;
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  assert(regions);
  assert(min_chunk_sz <= max_chunk_sz);

  Loops output;
  for (Loops::const_iterator i = regions->begin(); i != regions->end(); ++i) {
    utility::vector1<Loop> pieces;
    decompose(min_chunk_sz, max_chunk_sz, *i, &pieces);

    for (utility::vector1<Loop>::const_iterator j = pieces.begin(); j != pieces.end(); ++j)
      output.push_back(*j);
  }
  *regions = output;
}

void decompose(core::Size min_chunk_sz,
               core::Size max_chunk_sz,
               const protocols::loops::Loop& loop,
               utility::vector1<protocols::loops::Loop>* pieces) {
  using core::Size;
  using protocols::loops::Loop;
  using utility::vector1;
  assert(pieces);

  // Base case
  if (loop.length() <= max_chunk_sz) {
    pieces->push_back(loop);
    return;
  }

  // TODO(cmiles) check for off-by-1 errors
  Size midpoint = loop.start() + loop.length() / 2;
  Loop candidate_left(loop.start(), midpoint - 1);
  Loop candidate_right(midpoint, loop.stop());

  // Impossible to partition <loop> and satisfy min_chunk_sz requirement
  if (candidate_left.length() < min_chunk_sz ||
      candidate_right.length() < min_chunk_sz) {
    pieces->push_back(loop);
    return;
  }

  // TODO(cmiles) consider including additional information (e.g. secondary
  // structure, structural conservation) to inform pivot selection
  Size pivot = numeric::random::random_range(
      loop.start() + min_chunk_sz - 1,
      loop.stop() - min_chunk_sz + 1);

  // Recursively decompose the left- and right-hand sides of the pivot
  Loop left(loop.start(), pivot - 1);
  vector1<Loop> pieces_left;
  decompose(min_chunk_sz, max_chunk_sz, left, &pieces_left);

  Loop right(pivot, loop.stop());
  vector1<Loop> pieces_right;
  decompose(min_chunk_sz, max_chunk_sz, right, &pieces_right);

  // Update the output parameter
  std::copy(pieces_left.begin(), pieces_left.end(), std::back_inserter(*pieces));
  std::copy(pieces_right.begin(), pieces_right.end(), std::back_inserter(*pieces));
}

void find_regions_with_minimum_size(const core::sequence::SequenceAlignment& alignment,
                                    const core::Size unaligned_region_min_sz,
                                    protocols::loops::Loops* aligned_regions,
                                    protocols::loops::Loops* unaligned_regions) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;

  assert(aligned_regions);
  assert(unaligned_regions);

  Size pose_space_num_residues = 0;
  for (Size ii = 1; ii <= alignment.length(); ++ii)
    pose_space_num_residues = std::max(pose_space_num_residues, alignment.sequence(1)->resnum(ii));

  protocols::comparative_modeling::bounded_loops_from_alignment(
      pose_space_num_residues, unaligned_region_min_sz, alignment, unaligned_regions);

  *aligned_regions = unaligned_regions->invert(pose_space_num_residues);
}

void emit_intermediate(const core::pose::Pose& pose,
                       const std::string& filename) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  if (option[OptionKeys::abinitio::debug]())
    core::io::pdb::dump_pdb(pose, filename);
}

protocols::jd2::ThreadingJob const * const current_job() {
  using protocols::jd2::InnerJobCOP;
  using protocols::jd2::JobDistributor;
  using protocols::jd2::ThreadingJob;

  JobDistributor* jd2 = JobDistributor::get_instance();
  InnerJobCOP inner = jd2->current_job()->inner_job();
  return (ThreadingJob const * const) inner();
}

void add_cutpoint_variants(core::pose::Pose* pose) {
  const core::kinematics::FoldTree& tree(pose->fold_tree());
  for (core::Size i = 1; i <= pose->total_residue(); ++i) {
    if (!tree.is_cutpoint(i) || i >= (pose->total_residue() - 1))
      continue;

    core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::CUTPOINT_LOWER, i);
    core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::CUTPOINT_UPPER, i+1);
    TR.Debug << "Added cutpoint variants to residues " << i << ", " << i + 1 << std::endl;
  }
}

void remove_cutpoint_variants(core::pose::Pose* pose) {
  const core::kinematics::FoldTree& tree(pose->fold_tree());
  for (core::Size i = 1; i <= pose->total_residue(); ++i) {
    if (!tree.is_cutpoint(i) || i >= (pose->total_residue() - 1))
      continue;

    core::pose::remove_variant_type_from_pose_residue(*pose, core::chemical::CUTPOINT_LOWER, i);
    core::pose::remove_variant_type_from_pose_residue(*pose, core::chemical::CUTPOINT_UPPER, i+1);
    TR.Debug << "Removed cutpoint variants from residues " << i << ", " << i + 1 << std::endl;
  }
}

}  // namespace nonlocal
}  // namespace protocols
