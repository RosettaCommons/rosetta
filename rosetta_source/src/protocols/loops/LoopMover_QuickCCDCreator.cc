// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopMover_QuickCCDCreator.hh
/// @brief  Header for LoopMover_QuickCCDCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/loops/LoopMover_QuickCCDCreator.hh>

// Package Headers
#include <protocols/loops/LoopMover_QuickCCD.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace loops {

LoopMover_Perturb_QuickCCDCreator::~LoopMover_Perturb_QuickCCDCreator() {}

moves::MoverOP LoopMover_Perturb_QuickCCDCreator::create_mover() const {
  return new loops::LoopMover_Perturb_QuickCCD();
}

std::string LoopMover_Perturb_QuickCCDCreator::keyname() const {
  return "LoopMover_Perturb_QuickCCD";
}

} //namespace
} //namespace
