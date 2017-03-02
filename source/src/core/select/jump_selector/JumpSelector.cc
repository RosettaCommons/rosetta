// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpSelector.cc
/// @brief  The JumpSelector class identifies a subset of jumps from a Pose
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/jump_selector/JumpSelector.hh>

// Utility headers
#include <utility/tag/Tag.hh>

namespace core {
namespace select {
namespace jump_selector {

JumpSelector::JumpSelector() {}

JumpSelector::~JumpSelector() {}

/// @details Noop implementation in the base class in the case that a derived
/// class has no need to read data from an input tag
void JumpSelector::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap &
)
{}

} //namespace jump_selector
} //namespace select
} //namespace core
