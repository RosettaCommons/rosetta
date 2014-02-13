// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/DofUnlock.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/DofUnlock.hh>

// Package headers
#include <protocols/environment/ProtectedConformation.hh>

// Project headers
#include <core/pose/Pose.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.DofUnlock", basic::t_info);

namespace protocols {
namespace environment {

using core::environment::DofPassport;
using core::environment::DofPassportCOP;

  DofUnlock::DofUnlock( core::conformation::Conformation& conf, DofPassportCOP pass ) :
  conformation_( conf )
{
  conformation_.push_passport( pass );
}

DofUnlock::~DofUnlock(){
  conformation_.pop_passport();
}

} // environment
} // protocols
