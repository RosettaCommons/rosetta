// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/ExtendedPoseMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/simple_moves/ExtendedPoseMover.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Project headers
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/util.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

// Package headers
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>

// C/C++ headers
#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

// Option Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <utility/excn/Exceptions.hh>
#include <core/kinematics/Jump.hh>

namespace protocols {
namespace simple_moves {

ExtendedPoseMover::ExtendedPoseMover(const std::string& sequence,
                                     const std::string& residue_type_set)
    : sequence_(sequence), residue_type_set_(residue_type_set) {}

bool ExtendedPoseMover::valid() const {
  return sequence() != "";
}

void ExtendedPoseMover::apply(core::pose::Pose& pose) {
  // Ensure that this instance is in a valid state
  assert(valid());
  pose.clear();

  protocols::loops::Loops loops;
  core::pose::make_pose_from_sequence(pose, sequence_, residue_type_set_);
  protocols::loops::set_extended_torsions_and_idealize_loops(pose, loops);
}

std::string ExtendedPoseMover::get_name() const {
  return "ExtendedPoseMover";
}

const std::string& ExtendedPoseMover::sequence() const {
  return sequence_;
}

const std::string& ExtendedPoseMover::residue_type_set() const {
  return residue_type_set_;
}

void ExtendedPoseMover::sequence(const std::string& sequence) {
  sequence_ = sequence;
}

void ExtendedPoseMover::residue_type_set(const std::string& residue_type_set) {
  residue_type_set_ = residue_type_set;
}

protocols::moves::MoverOP ExtendedPoseMover::clone() const {
  return new protocols::simple_moves::ExtendedPoseMover(*this);
}

protocols::moves::MoverOP ExtendedPoseMover::fresh_instance() const {
  return new protocols::simple_moves::ExtendedPoseMover();
}

void ExtendedPoseMover::parse_my_tag(const utility::tag::TagCOP tag,
	                                         basic::datacache::DataMap&,
                                     const protocols::filters::Filters_map&,
                                     const protocols::moves::Movers_map&,
                                     const core::pose::Pose&) {
  // required options
  if (tag->hasOption("sequence")) sequence(tag->getOption<string>("sequence")); 
  else if (tag->getOption<bool>("use_fasta", false )){
    sequence( core::sequence::read_fasta_file_return_str(
        basic::options::option[ basic::options::OptionKeys::in::file::fasta ]()[1] ) );
  }
  else throw utility::excn::EXCN_RosettaScriptsOption("Failed to specify required option `sequence` or fasta file");

  // additional options
  if (tag->hasOption("residue_type_set"))
    residue_type_set(tag->getOption<string>("residue_type_set"));
}

}  // namespace simple_moves
}  // namespace protocols
