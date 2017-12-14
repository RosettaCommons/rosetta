// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/execute.hh
/// @brief  Various functions to run external executables
/// @author Sergey Lyskov
///
/// Note: this code is placed in to 'basic' level instead of utility because we want to have access to Tracers so verbosity of all external calls could be controlled


#include <basic/execute.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>


#include <iostream>

namespace basic {

static basic::Tracer TR("basic.execute");

#ifdef WIN32
#define popen  _popen
#define pclose _pclose
#endif

ExecutionResult execute(std::string const & message, std::string const & command_line, bool terminate_on_failure, bool silent)
{
	ExecutionResult r;

	if ( !silent )  TR << message << ' ' << command_line << std::endl;  // TR.flush();

	FILE *pipe = popen( (command_line+ " 2>&1" ).c_str(), "r"); // we need to redirect std error in to std output so it will got captured by pipe...

	if ( pipe ) {
		char buffer[256];

		while ( !feof(pipe) ) {
			if ( fgets(buffer, 256, pipe) != nullptr ) r.output += buffer;
		}

		r.result = pclose(pipe);
	} else {
		r.result = 1;
		r.output = "Could not open pipe!";
	}

	if ( terminate_on_failure && r.result ) utility_exit_with_message("basic.execute encounter error while runnning " + command_line+ '\n'+ r.output +"\nExiting...");
	return r;
}


} // namespace basic
