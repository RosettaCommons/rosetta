// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/EnvBeginMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/EnvBeginMover.hh>

// Package headers
#include <protocols/environment/Environment.hh>

// Project headers
#include <core/pose/Pose.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.movers.EnvBeginMover", basic::t_info);

namespace protocols {
namespace environment {

EnvBeginMover::EnvBeginMover():
  Mover()
{}

EnvBeginMover::EnvBeginMover( EnvironmentOP env ):
  env_( env )
{}

EnvBeginMover::~EnvBeginMover() {}

void EnvBeginMover::apply( Pose& pose ) {
  pose = env_->start( pose );
}

std::string EnvBeginMover::get_name() const {
  return "EnvBeginMover("+env_->name()+")";
}

moves::MoverOP EnvBeginMover::clone() const {
  return new EnvBeginMover( *this );
}



} // environment
} // protocols
