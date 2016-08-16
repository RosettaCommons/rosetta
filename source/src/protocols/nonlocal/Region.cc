// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/Region.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/Region.hh>

// Project headers
#include <core/types.hh>

namespace protocols {
namespace nonlocal {

/// @details Auto-generated virtual destructor
Region::~Region() {}

Region::Region(core::Size start_pos, core::Size stop_pos)
: start_(start_pos), stop_(stop_pos) {}

core::Size Region::start() const {
	return start_;
}

core::Size Region::stop() const {
	return stop_;
}

core::Size Region::length() const {
	return increasing() ? stop() - start() + 1 : start() - stop() + 1;
}

bool Region::increasing() const {
	return start() <= stop();
}

bool Region::decreasing() const {
	return stop() <= start();
}

}  // namespace nonlocal
}  // namespace protocols
