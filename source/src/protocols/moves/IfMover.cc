// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/IfMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/moves/IfMover.hh>
#include <protocols/moves/IfMoverCreator.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {

static thread_local basic::Tracer TR( "protocols.moves.IfMover" );

std::string IfMoverCreator::mover_name() {
  return "If";
}

std::string IfMoverCreator::keyname() const {
  return mover_name();
}

protocols::moves::MoverOP IfMoverCreator::create_mover() const {
  return new IfMover();
}

void IfMover::apply(core::pose::Pose& pose) {
  moves::MoverStatus status;

  if( filter_->apply(pose) ){
    true_mover_->apply(pose);
    status = true_mover_->get_last_move_status();
  } else {
    false_mover_->apply(pose);
    status = false_mover_->get_last_move_status();
  }

  // update this mover's status
  protocols::moves::Mover::set_last_move_status(status);
}

core::pose::PoseOP IfMover::get_additional_output_true_mover() {
  return true_mover_->get_additional_output();
}

core::pose::PoseOP IfMover::get_additional_output_false_mover() {
  return false_mover_->get_additional_output();
}

// backwards compatibility
core::pose::PoseOP IfMover::get_additional_output() {
  return get_additional_output_true_mover();
}

std::string IfMover::get_name() const {
  return IfMoverCreator::mover_name();
}

void IfMover::parse_my_tag( utility::tag::TagCOP tag,
                            basic::datacache::DataMap &,
                            protocols::filters::Filters_map const &filters,
                            protocols::moves::Movers_map const &movers,
                            core::pose::Pose const & ) {
  using namespace protocols::filters;

  TR<<"If mover\n";
  std::string const true_mover_name( tag->getOption< std::string >( "true_mover_name" ));
  std::string const false_mover_name( tag->getOption< std::string >( "false_mover_name", "null" ));
  std::string const filter_name( tag->getOption< std::string >( "filter_name" ) );

	/// see: protocols/moves/util.hh
  filter_ = find_filter_or_die(filter_name, tag, filters);
  true_mover_ = find_mover_or_die(true_mover_name, tag, movers);
  false_mover_ = find_mover_or_die(false_mover_name, tag, movers);

  TR << "with true_mover \"" << true_mover_name
     << "\" and false_mover \"" << false_mover_name
     << "\" filter \"" << filter_name
     << std::endl;
}

} //moves
} //protocols

