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
#include <cmath>
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/util.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/nonlocal/CutFinder.hh>
#include <protocols/nonlocal/NLGrouping.hh>

namespace protocols {
namespace nonlocal {

static basic::Tracer TR("protocols.nonlocal.util");
static numeric::random::RandomGenerator RG(156144120);

void nonlocal_groupings_from_alignment(utility::vector1<NLGrouping>* groupings) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;
  using core::id::SequenceMapping;
  using core::pose::PoseOP;
  using core::sequence::SequenceAlignment;
  using core::sequence::SequenceAlignmentOP;
  using protocols::loops::Loops;
  using std::string;
  using utility::vector1;

  // Assumes multiple alignments (if present) are contained in a single file
  string alignment_file = option[in::file::alignment]()[1];
  vector1<SequenceAlignment> alignments = core::sequence::read_aln(option[cm::aln_format](), alignment_file);
  assert(alignments.size() > 0);

  SequenceAlignmentOP alignment = alignments[1].clone();
  PoseOP template_pose = core::import_pose::pose_from_pdb(option[in::file::template_pdb]()[1]);

  // Identify stretches of aligned and unaligned residues in the alignment.
  // Limit the length of the aligned regions to enhance conformational sampling.
  FragmentIO io;
  FragSetOP fragments = io.read_data(option[in::file::frag3]());
  core::Size min_chunk_sz = fragments->max_frag_length();
  core::Size max_chunk_sz = option[OptionKeys::nonlocal::max_chunk_size]();

  Loops aligned_regions, unaligned_regions;
  find_regions(*alignment, min_chunk_sz, &aligned_regions, &unaligned_regions);
  limit_chunk_size(min_chunk_sz, max_chunk_sz, &aligned_regions);

  TR.Debug << "Aligned:" << std::endl << aligned_regions << std::endl;
  TR.Debug << "Unaligned:" << std::endl << unaligned_regions << std::endl;

  // Generate non-local groupings by iterating over the regions identified
  // above. Backbone torsions and coordinates are taken from the template
  NLGrouping grouping;
  generate_nonlocal_grouping(aligned_regions,
                             unaligned_regions,
                             alignment->sequence_mapping(1, 2),
                             *template_pose,
                             &grouping);

  // update the output parameter
  groupings->push_back(grouping);
  TR << grouping.provenance() << std::endl;
}

/// @brief Enforces restrictions on min/max chunk size in <regions> by applying
/// recursive decomposition
void limit_chunk_size(core::Size min_chunk_sz,
                      core::Size max_chunk_sz,
                      protocols::loops::Loops* regions) {
  using core::Size;
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  assert(regions);
  assert(min_chunk_sz > 0);
  assert(max_chunk_sz > 0);
  assert(min_chunk_sz <= max_chunk_sz);

  Loops output;
  for (Loops::const_iterator i = regions->begin(); i != regions->end(); ++i) {
    const Loop& loop = *i;

    // TODO(cmiles) trouble if we encounter aligned regions w/ length < min_chunk_sz.
    // Gap sampling extension partially mitigates unaligned case.
    if (loop.length() <= max_chunk_sz) {
      output.push_back(loop);
      continue;
    }

    // Recursively decompose <loop> until each piece has length <= max_length
    utility::vector1<Loop> pieces;
    decompose(min_chunk_sz, max_chunk_sz, loop, &pieces);

    for (utility::vector1<Loop>::const_iterator j = pieces.begin(); j != pieces.end(); ++j) {
      const Loop& loop = *j;
      if (loop.length() >= min_chunk_sz)
        output.push_back(loop);
    }
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
  Loop candidate_left(loop.start(), midpoint);
  Loop candidate_right(midpoint, loop.stop());

  // Impossible to satisfy min_chunk_sz requirement
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
  Loop left(loop.start(), pivot);
  vector1<Loop> pieces_left;
  decompose(min_chunk_sz, max_chunk_sz, left, &pieces_left);

  Loop right(pivot + 1, loop.stop());
  vector1<Loop> pieces_right;
  decompose(min_chunk_sz, max_chunk_sz, right, &pieces_right);

  // Update the output parameter
  std::copy(pieces_left.begin(), pieces_left.end(), std::back_inserter(*pieces));
  std::copy(pieces_right.begin(), pieces_right.end(), std::back_inserter(*pieces));
}

void find_regions(const core::sequence::SequenceAlignment& alignment,
                  const core::Size unaligned_region_min_sz,
                  protocols::loops::Loops* aligned_regions,
                  protocols::loops::Loops* unaligned_regions) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;

  assert(aligned_regions);
  assert(unaligned_regions);
  assert(unaligned_region_min_sz > 0);

  Size pose_space_num_residues = 0;
  for (Size ii = 1; ii <= alignment.length(); ++ii)
    pose_space_num_residues = std::max(pose_space_num_residues, alignment.sequence(1)->resnum(ii));

  *unaligned_regions = protocols::comparative_modeling::loops_from_alignment( pose_space_num_residues, alignment, unaligned_region_min_sz);
  *aligned_regions   = unaligned_regions->invert(pose_space_num_residues);
}

void generate_nonlocal_grouping(const protocols::loops::Loops& aligned_regions,
                                const protocols::loops::Loops& unaligned_regions,
                                const core::id::SequenceMapping& mapping,
                                const core::pose::Pose& template_pose,
                                NLGrouping* grouping) {
  using core::Real;
  using core::Size;
  assert(grouping);

  bool unalignedBegin = false, unalignedEnd = false;
  if (unaligned_regions.size() > 0) {
    unalignedBegin = unaligned_regions[1].start() < aligned_regions[1].start();
    unalignedEnd = unaligned_regions[unaligned_regions.size()].start() > aligned_regions[aligned_regions.size()].start();
  }
  TR.Debug << "Unaligned begin: " << unalignedBegin << std::endl;
  TR.Debug << "Unaligned end: " << unalignedEnd << std::endl;

  for (Size i = 1; i <= aligned_regions.num_loop(); ++i) {
    NLFragmentGroup group;

    const protocols::loops::Loop& region = aligned_regions[i];
    for (Size j = region.start(); j <= region.stop(); ++j) {
      Size mapped_index = mapping[j];

      // torsions
      Real phi = template_pose.phi(mapped_index);
      Real psi = template_pose.psi(mapped_index);
      Real omega = template_pose.omega(mapped_index);

      // CA coordinates
      const core::PointPosition& coords = template_pose.xyz(core::id::NamedAtomID("CA", mapped_index));
      group.add_entry(NLFragment(j, phi, psi, omega, coords.x(), coords.y(), coords.z()));
    }
    grouping->add_group(group);
  }
  grouping->sort();
}

bool has_chainbreaks(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  protocols::loops::Loops loops =
      protocols::comparative_modeling::pick_loops_chainbreak(pose, option[cm::min_loop_size]());
  loops.verify_against(pose);
  return loops.num_loop() > 0;
}

void emit_intermediate(const core::pose::Pose& pose,
                       const std::string& filename) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  if (option[OptionKeys::abinitio::debug]())
    core::io::pdb::dump_pdb(pose, filename);
}

}  // namespace nonlocal
}  // namespace protocols
