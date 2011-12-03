// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopMover_KICCreator.hh
/// @brief  Header for LoopMover_KICCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/loops/LoopMover_KICCreator.hh>

// Package Headers
#include <protocols/loops/LoopMover_KIC.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace loops {

LoopMover_Perturb_KICCreator::~LoopMover_Perturb_KICCreator() {}

moves::MoverOP LoopMover_Perturb_KICCreator::create_mover() const {
  return new loops::LoopMover_Perturb_KIC();
}

std::string LoopMover_Perturb_KICCreator::keyname() const {
  return "LoopMover_Perturb_KIC";
}

} //namespace
} //namespace
