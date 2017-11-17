// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/environment/EnvExcn.hh
/// @brief A special kind of exception for the bad stuff that can happen in exceptions.
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
#include <sstream>

// ObjexxFCL Headers

namespace core {
namespace environment {

class EXCN_Env_Exception : public utility::excn::Exception {
	typedef utility::excn::Exception Parent;
public:
	EXCN_Env_Exception(char const *file, int line, EnvCoreCAP env_ap ) : Parent(file, line, "") {
		EnvCoreCOP env = env_ap.lock(); // avoids throwing weak_ptr exception
		std::ostringstream elab_msg;

		if ( !env ) {
			elab_msg << "Environment exception in unknown environment!";

		} else {
			elab_msg << "Environment exception in environment: '" << env->name();

			EnvCoreCAP superenv = env->superenv();
			while ( ! superenv.expired() ) {
				EnvCoreCOP superenv_op = superenv.lock();
				if ( !superenv_op ) break;
				elab_msg << superenv_op->name() << ">";
				superenv = superenv_op->superenv();
			}

			elab_msg << "'";
		}

		add_msg( elab_msg.str() );
	}

protected:
	//EXCN_Env_Exception(char const *file, int line) : Parent(file, line, "") {};
};

} // environment
} // core

#endif //INCLUDED_core_environment_EnvExcnBase_hh
