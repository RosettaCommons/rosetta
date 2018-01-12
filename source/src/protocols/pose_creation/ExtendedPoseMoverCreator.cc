// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/ExtendedPoseMoverCreator.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/pose_creation/ExtendedPoseMoverCreator.hh>

// C/C++ headers
#include <string>

// Package headers
#include <protocols/pose_creation/ExtendedPoseMover.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace pose_creation {

// XRW TEMP protocols::moves::MoverOP ExtendedPoseMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new protocols::pose_creation::ExtendedPoseMover() );
// XRW TEMP }

// XRW TEMP std::string ExtendedPoseMoverCreator::keyname() const {
// XRW TEMP  return mover_name();
// XRW TEMP }

// XRW TEMP std::string ExtendedPoseMoverCreator::mover_name() {
// XRW TEMP  return "ExtendedPoseMover";
// XRW TEMP }

}  // namespace pose_creation
}  // namespace protocols
