// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   init_util.hh
/// @brief  initialization subroutines for Benchmark
/// @author

#ifndef INCLUDED_apps_benchmark_init_util_hh
#define INCLUDED_apps_benchmark_init_util_hh

#include <devel/init.hh>

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// @brief For unit tests only.  Constructs a mock argc and argv given a stl string.
inline
char**
create_pseudo_commandline( std::string const & cmdin, int & pseudo_argc )
{
	using namespace std;
	pseudo_argc = 0;

	vector< int > string_start;
	vector< int > string_length;
	bool in_the_middle_of_a_string = false;
	for ( unsigned int ii = 0; ii < cmdin.size(); ++ii ) {
		if ( cmdin[ ii ] != ' ' ) {
			if ( in_the_middle_of_a_string ) {
				++string_length[ pseudo_argc - 1 ];
			} else {
				string_start.push_back( ii );
				string_length.push_back( 1 );
				in_the_middle_of_a_string = true;
				++pseudo_argc;
			}
		} else {
			in_the_middle_of_a_string = false;
		}
	}

	vector< string > cmdline;
	cmdline.reserve( pseudo_argc );
	for ( int ii = 0; ii < pseudo_argc; ++ii ) {
		cmdline.push_back( cmdin.substr( string_start[ ii ], string_length[ ii ] ));
	}

	char** pseudo_argv = new char * [ pseudo_argc ];
	for ( int ii = 0; ii < pseudo_argc; ++ii ) {
		pseudo_argv[ ii ] = new char[ cmdline[ ii ].size() + 1 ];
		strncpy(pseudo_argv[ii], cmdline[ii].c_str(), cmdline[ii].size());
		pseudo_argv[ii][ cmdline[ii].size() ] = '\0';
	}

	return pseudo_argv;
}

// @brief For unit tests only.  Deallocates the pseudo command line created by a call to
// create_pseudo_commandline
inline
void
destroy_pseudo_commandline( int pseudo_argc, char** pseudo_argv )
{
	for ( int ii = 0; ii < pseudo_argc; ++ii ) {
		delete [] pseudo_argv[ ii ];
	}
	delete [] pseudo_argv;
}


/// @brief For unit tests only.  Creates an argc/argv pair, calls init() and deletes argv.
/// can be used to init option system more that once, providing user ability to change
/// options on the fly.
inline void core_init_from_string( std::string const & commandline )
{
	extern char ** command_line_argv;  /// We need original filename to be first argument

	int pseudo_argc;
	char** pseudo_argv = create_pseudo_commandline( std::string(command_line_argv[0]) + " "
		+ commandline, pseudo_argc );
	devel::init( pseudo_argc, pseudo_argv );
	destroy_pseudo_commandline( pseudo_argc, pseudo_argv );
}

/// @brief For unit tests only. Re-init option system.
/// Command line will be = old command line + function argument.
inline void core_init_with_additional_options( std::string const & commandline_in )
{
	extern int command_line_argc; extern char ** command_line_argv;

	std::string commandline(" ");
	for ( int i=1; i<command_line_argc; i++ ) {
		commandline = commandline + command_line_argv[i] + " ";
	}
	commandline = commandline + commandline_in;
	//std::cout << "core_init_with_additional_options: " << commandline.c_str() << "\n";

	core_init_from_string( commandline );
}


#endif
