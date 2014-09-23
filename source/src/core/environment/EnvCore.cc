// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/environment/EnvCore.cc
/// @author Justin Porter

// Unit Headers
#include <core/environment/EnvCore.hh>

// Package headers
#include <core/environment/DofPassport.hh>

// Project headers
#include <utility/excn/Exceptions.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static thread_local basic::Tracer tr( "protocols.environment.EnvCore", basic::t_info );

namespace core {
namespace environment {

core::Size EnvCore::current_maximum_id_ = 0;

EnvCore::EnvCore( std::string const& name ):
  name_( name ),
  id_( generate_id() )
{}

EnvCore::~EnvCore() {}

std::string const& EnvCore::name() const{
  return name_;
}

DofPassportOP EnvCore::issue_passport( std::string const& mover_name ) const{
	// This method doesn't register the passport with the mover because that would
	// require including movers, which violates the inclusion hierarchy.
  return DofPassportOP( new DofPassport( mover_name, id() ) );
}

EnvCoreCAP EnvCore::superenv() const{
  return superenv_;
}

void EnvCore::set_superenv( EnvCoreCAP env ){
  if( superenv_.expired() || utility::pointer::equal(superenv_, env) ){
    superenv_ = env;
  } else {
    throw utility::excn::EXCN_Msg_Exception( "Superenvironment already assigned.");
  }
}

core::Size const& EnvCore::id() const {
  return id_;
}

core::Size EnvCore::generate_id(){
  current_maximum_id_ += 1;
  return current_maximum_id_;
}


} // environment
} // core
