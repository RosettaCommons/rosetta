// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/util.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_UTIL_HH_
#define PROTOCOLS_NONLOCAL_UTIL_HH_

// C/C++ headers
#include <string>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/nonlocal/NLGrouping.hh>

namespace protocols {
namespace nonlocal {

/// @brief Creates a set of non-local groupings from a threading model.
/// The sequence alignment and template pdb are read from the filenames
/// specified in the options -in:file:alignment and -in:file:template_pdb,
/// respectively.
void nonlocal_groupings_from_alignment(utility::vector1<NLGrouping>* groupings);

// -- Utility methods -- not to be called directly
void find_regions(const core::sequence::SequenceAlignment& alignment,
                  const core::Size num_residues,
                  protocols::loops::Loops* aligned_regions,
                  protocols::loops::Loops* unaligned_regions);

void generate_nonlocal_grouping(const protocols::loops::Loops& aligned_regions,
                                const protocols::loops::Loops& unaligned_regions,
                                const core::id::SequenceMapping& mapping,
                                const core::pose::Pose& template_pose,
                                NLGrouping* grouping);

void read_from_template(const protocols::loops::Loops& regions,
                        const core::id::SequenceMapping& mapping,
                        const core::pose::Pose& template_pose,
                        NLGrouping* grouping);

/// @brief Examines the contents of <regions>. If the length of any element
/// exceeds -nonlocal:max_chunk_size, recursively decomposes that region into
/// a series of pieces each less than the threshold.
void limit_chunk_size(core::Size min_chunk_sz,
                      core::Size max_chunk_sz,
                      protocols::loops::Loops* regions);

/// @brief Recursively decomposes <loop> into a series of <pieces>, each having
/// length less than or equal to <max_length>.
void decompose(core::Size min_chunk_sz,
               core::Size max_chunk_sz,
               const protocols::loops::Loop& loop,
               utility::vector1<protocols::loops::Loop>* pieces);

/// @brief Returns true if <pose> has chainbreaks, false otherwise.
/// Non-const because protocols/comparative_modeling/util.cc's methods are
/// unfortunately non-const.
bool has_chainbreaks(core::pose::Pose& pose);

/// @brief If -abinitio:debug is enabled, writes <pose> to <file>.
void emit_intermediate(const core::pose::Pose& pose, const std::string& file);

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_UTIL_HH_
