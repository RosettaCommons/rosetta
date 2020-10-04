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


#include <utility/string_util.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>

#include <iostream>

#include <sys/wait.h>

namespace basic {

static basic::Tracer TR("basic.execute");

#ifdef WIN32
#define popen  _popen
#define pclose _pclose
#endif

ExecutionResult execute(std::string const & message, std::string const & command, std::vector<std::string> const & args, bool terminate_on_failure, bool silent)
{
	if ( !silent )  TR << message << ' ' << command << ' ' << utility::to_string(args) << std::endl;  // TR.flush();

	ExecutionResult r {1, ""};

	int pipe_fd[2];

	if ( pipe(pipe_fd) == -1 ) utility_exit_with_message("Error creating pipe...");

	pid_t pid = fork();

	if ( pid ) { // parent
		close(pipe_fd[1]); // parent does not write

		int const buffer_size = 256;
		char buffer[buffer_size];

		while ( int count = read(pipe_fd[0], &buffer, sizeof(buffer)) ) r.output.append(buffer, count);

		waitpid(-1, &r.result, 0);

		close(pipe_fd[0]);
	} else { // child
		close(pipe_fd[0]); // child does not read

		dup2(pipe_fd[1], 1); // redirecting stdout
		dup2(pipe_fd[1], 2); // redirecting stderr

		std::vector<char const *> c_args;
		c_args.push_back(command.c_str()); // first argument should be program name
		for ( auto const & s : args ) c_args.push_back(s.c_str());
		c_args.push_back(nullptr);

		//execvp(command.c_str(), static_cast<char * const *>( c_args.data() ) );
		execvp(command.c_str(), (char * const *)c_args.data() );

		utility_exit_with_message("ERROR: execvp failed, command: '" + command + "' args: " + utility::to_string(args) );
	}

	if ( terminate_on_failure && r.result ) utility_exit_with_message("basic.execute encounter error while runnning " + command + ' ' + utility::to_string(args) + '\n'+ r.output +"\nExiting...");
	return r;
}


} // namespace basic
