// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/EnvEndMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/EnvEndMover.hh>

// Package headers
#include <protocols/environment/Environment.hh>

// Project headers
#include <core/pose/Pose.hh>


// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.movers.EnvEndMover", basic::t_info);

namespace protocols {
namespace environment {

using namespace core::environment;

EnvEndMover::EnvEndMover( EnvironmentOP env ):
  env_( env )
{}

EnvEndMover::~EnvEndMover() {}

void EnvEndMover::apply( core::pose::Pose& pose ) {
  pose = env_->end( pose );
}

std::string EnvEndMover::get_name() const {
  return "EnvEndMover("+env_->name()+")";
}

moves::MoverOP EnvEndMover::clone() const {
  return new EnvEndMover( *this );
}

} // environment
} // protocols
