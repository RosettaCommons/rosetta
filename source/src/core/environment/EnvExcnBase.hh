// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/environment/EnvExcn.hh
/// @brief A class that is used to express a mover-specific DoF-unlock on a ProtectedPose. Its destruction expresses a re-locking.
///
/// @author Justin Porter

#ifndef INCLUDED_core_environment_EnvExcnBase_hh
#define INCLUDED_core_environment_EnvExcnBase_hh

// Unit Headers

// Package headers
#include <core/environment/EnvCore.hh>

#include <utility/excn/Exceptions.hh>

// Project headers

// C++ Headers
#include <string>

// ObjexxFCL Headers

namespace core {
namespace environment {

class EXCN_Env_Exception : public utility::excn::EXCN_Msg_Exception {
	typedef utility::excn::EXCN_Msg_Exception Parent;
public:
	EXCN_Env_Exception( std::string const& msg, EnvCoreCAP env ) {
		std::ostringstream elab_msg;
    elab_msg << msg << " Logged by Environment: '";

    EnvCoreCAP superenv = env->superenv().get();
    while( superenv != 0 ){
      elab_msg << superenv->name() << ">";
      superenv = superenv->superenv().get();
    }

		Parent( elab_msg.str() );
	}
protected:
  EXCN_Env_Exception() : Parent() {};
};

} // environment
} // core

#endif //INCLUDED_core_environment_EnvExcnBase_hh
