// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/ExtendedPoseMoverCreator.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/ExtendedPoseMoverCreator.hh>

// C/C++ headers
#include <string>

// Project headers
#include <protocols/moves/Mover.fwd.hh>

// Package headers
#include <protocols/nonlocal/ExtendedPoseMover.hh>

namespace protocols {
namespace nonlocal {

protocols::moves::MoverOP ExtendedPoseMoverCreator::create_mover() const {
  return new ExtendedPoseMover();
}

std::string ExtendedPoseMoverCreator::keyname() const {
  return mover_name();
}

std::string ExtendedPoseMoverCreator::mover_name() {
  return "ExtendedPoseMover";
}

}  // namespace nonlocal
}  // namespace protocols
