// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/init_util.hh
/// @brief  initialization subroutines for unit testing
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_protocols_init_util_HH
#define INCLUDED_protocols_init_util_HH

#include <core/init_util.hh>

#include <protocols/init/init.hh>

#include <string>
#include <fstream>
#include <iostream>
#include <vector>


/// @brief For unit tests only.  Creates an argc/argv pair, calls init() and deletes argv.
/// can be used to init option system / factories more that once, providing user ability to change
/// options on the fly while guaranteeing factory initialization
inline void protocols_init_from_string( std::string const & commandline )
{
	extern char ** command_line_argv;  /// We need original filename to be first argument

	int pseudo_argc;
	char** pseudo_argv = create_pseudo_commandline(
		std::string(command_line_argv[0]) + " "
		+ commandline, pseudo_argc );
	protocols::init::init( pseudo_argc, pseudo_argv );
	destroy_pseudo_commandline( pseudo_argc, pseudo_argv );
	initialize_rng();
}

/// @brief For unit tests only. Re-init option system.
/// Command line will be = old command line + function argument.
inline void protocols_init_with_additional_options( std::string const & commandline_in )
{
	protocols_init_from_string( get_commandline_with_additional_options(commandline_in) );
}


// @brief Analog of protocols::init::init() for unit test suite.
//
inline void protocols_init()
{
	extern int command_line_argc;
	extern char ** command_line_argv;

	if( command_line_argc > 1 ) {
		protocols::init::init(command_line_argc, command_line_argv);
		initialize_rng();
	}
	else {
		std::string commandline = "core.test -mute all";
		std::string db_cmdline = append_db_to_commandline( commandline );
		protocols_init_from_string( db_cmdline );
	}
}

#endif
