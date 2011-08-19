// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/ExtendedPoseMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/ExtendedPoseMover.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Project headers
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>

// C/C++ headers
#include <cassert>
#include <string>

namespace protocols {
namespace nonlocal {

ExtendedPoseMover::ExtendedPoseMover(const std::string& sequence,
                                     const std::string& residue_type_set)
    : sequence_(sequence), residue_type_set_(residue_type_set) {}

bool ExtendedPoseMover::valid() const {
  return sequence() != "";
}

// -- Mover -- //
void ExtendedPoseMover::apply(core::pose::Pose& pose) {
  // Ensure that this instance is in a valid state
  assert(valid());

  // Completely wipe out any existing contents in <pose>
  pose.clear();

  protocols::loops::Loops loops;
  core::pose::make_pose_from_sequence(pose, sequence_, residue_type_set_);
  protocols::loops::set_extended_torsions_and_idealize_loops(pose, loops);
}

std::string ExtendedPoseMover::get_name() const {
  return "ExtendedPoseMover";
}

// -- Accessors -- //
const std::string& ExtendedPoseMover::sequence() const {
  return sequence_;
}

const std::string& ExtendedPoseMover::residue_type_set() const {
  return residue_type_set_;
}

// -- Mutators -- //
void ExtendedPoseMover::sequence(const std::string& sequence) {
  sequence_ = sequence;
}

void ExtendedPoseMover::residue_type_set(const std::string& residue_type_set) {
  residue_type_set_ = residue_type_set;
}


// -- RosettaScripts -- //
protocols::moves::MoverOP ExtendedPoseMover::clone() const {
  return new ExtendedPoseMover(*this);
}

protocols::moves::MoverOP ExtendedPoseMover::fresh_instance() const {
  return new ExtendedPoseMover();
}

void ExtendedPoseMover::parse_my_tag(const utility::tag::TagPtr tag,
                                     protocols::moves::DataMap&,
                                     const protocols::filters::Filters_map&,
                                     const protocols::moves::Movers_map&,
                                     const core::pose::Pose&) {
  // required options
  if (!tag->hasOption("sequence"))
    utility_exit_with_message("Failed to specify required option `sequence`");

  sequence(tag->getOption<string>("sequence"));

  // additional options
  if (tag->hasOption("residue_type_set"))
    residue_type_set(tag->getOption<string>("residue_type_set"));
}

}  // namespace nonlocal
}  // namespace protocols
