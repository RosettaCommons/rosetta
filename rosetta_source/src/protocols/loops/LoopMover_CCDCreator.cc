// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopMover_CCDCreator.hh
/// @brief  Header for LoopMover_CCDCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/loops/LoopMover_CCDCreator.hh>

// Package Headers
#include <protocols/loops/LoopMover_CCD.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace loops {

LoopMover_Perturb_CCDCreator::~LoopMover_Perturb_CCDCreator() {}

moves::MoverOP LoopMover_Perturb_CCDCreator::create_mover() const {
  return new loops::LoopMover_Perturb_CCD();
}

std::string LoopMover_Perturb_CCDCreator::keyname() const {
  return "LoopMover_Perturb_CCD";
}




LoopMover_Refine_CCDCreator::~LoopMover_Refine_CCDCreator() {}

moves::MoverOP LoopMover_Refine_CCDCreator::create_mover() const {
  return new loops::LoopMover_Refine_CCD();
}

std::string LoopMover_Refine_CCDCreator::keyname() const {
  return "LoopMover_Refine_CCD";
}



} //namespace
} //namespace
