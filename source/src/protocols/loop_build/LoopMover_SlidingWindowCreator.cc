// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loop_build/LoopMover_SlidingWindowCreator.hh
/// @brief  Header for LoopMover_SlidingWindowCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/loop_build/LoopMover_SlidingWindowCreator.hh>

// Package Headers
#include <protocols/loop_build/LoopMover_SlidingWindow.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace loop_build {

LoopMover_SlidingWindowCreator::~LoopMover_SlidingWindowCreator() {}

moves::MoverOP LoopMover_SlidingWindowCreator::create_mover() const {
	return moves::MoverOP( new loop_build::LoopMover_SlidingWindow() );
}

std::string LoopMover_SlidingWindowCreator::keyname() const {
	return "LoopMover_SlidingWindow";
}

} //namespace
} //namespace
