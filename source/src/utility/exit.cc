// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/exit.cc
/// @brief  Program exit functions and macros
/// @author David Kim (dekim@u.washington.edu)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)

#ifdef USEMPI
#include <mpi.h>
#endif
// Unit headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/backtrace.hh>
#include <utility/CSI_Sequence.hh>

// C++ headers
#include <cassert>
#include <cstdlib>
#include <iostream>

// Boinc headers
#ifdef BOINC
#include <utility/io/izstream.hh>
#include <boinc_api.h>
#endif

#include <string>
#include <vector>
#include <sstream>

#ifndef _WIN32
#include <unistd.h>
#endif

#ifdef _WIN32
#include <io.h>
//__declspec(dllexport) int isatty(int fd) { return _isatty(fd); }
#endif

namespace utility {


struct UtilityExitException : excn::Exception //noboday should be allowed to throw this... that's why its privately hidden in this modul...
{
	UtilityExitException(char const * file, int line, std::string const& msg) : Exception(file, line, "[ ERROR ] UtilityExitException\nERROR: " + msg + "\n\n")
	{}
};



/// Place holder for 'end-action' of utility::exit(…)
static void (*main_exit_callback)(void) = nullptr;

void set_main_exit_callback( UtilityExitCallBack my_callback )
{
	main_exit_callback = my_callback;
}

/// Array to hold all additional exit-callbacks
std::vector< UtilityExitCallBack > & get_all_exit_callbacks()
{
	static auto * all_CB = new std::vector< UtilityExitCallBack >;
	return *all_CB;
}

void add_exit_callback( UtilityExitCallBack cb)
{
	get_all_exit_callbacks().push_back( cb );
}

void remove_exit_callback( UtilityExitCallBack cb )
{
	for ( auto it=get_all_exit_callbacks().begin(); it < get_all_exit_callbacks().end(); ++it ) {
		if ( (*it) == cb ) {
			get_all_exit_callbacks().erase(it);
			break;
		}
	}
}



/// @brief Exit with file + line + message + optional status
void
exit(
	char const * file,
	int const line,
	std::string const & message,
	int const status
)
{

	// Calling all preset exit-callback's
	for ( auto it=get_all_exit_callbacks().begin(); it < get_all_exit_callbacks().end(); ++it ) {
		(*it)();
	}

	std::ostringstream oss;
	if ( ! message.empty() ) oss << "\n" << "ERROR: " << message << "\n";
	oss << "ERROR:: Exit from: " << file << " line: " << line << "\n";
	std::string failure_message = oss.str();
	maybe_throw_on_next_assertion_failure( failure_message.c_str() );

	std::cerr << CSI_Reset() << CSI_Red() << CSI_Bold();
	std::cerr << failure_message << std::flush;
	std::cerr << CSI_Reset();
	print_backtrace( message.c_str() );
	std::cerr.flush();

#ifndef BOINC
	// why did this get placed here basically skipping the logic below?!
	throw UtilityExitException(file, line,  message);
#endif

#ifdef USEMPI
	MPI_Abort( MPI_COMM_WORLD, 911 );
#endif

#ifdef BOINC

	// check if there are results, if so, return success
	bool hasresults = false;
	// Quick hack, since we can't access the options assume
	// the result file will always be named default.out and thaat
	// it should be gzipped for now.
	utility::io::izstream data( "default.out" );
	if ( !data ){
		std::cerr << "BOINC:: Error reading and gzipping output datafile: default.out" << std::endl; std::cerr.flush();
		boinc_finish( status );
	}
	std::string tmpline;
	getline( data, tmpline ); // sequence line
	getline( data, tmpline ); // score line
	while( getline(data,tmpline) ) {
		if ( tmpline.substr(0,7) == "SCORE: " ) {
			hasresults = true;
			break;
		}
	}
	data.close();
	if (hasresults) {
		utility::file::gzip( "default.out", true );
		boinc_finish( 0 );
	} else {
		boinc_finish( status );
	}

#else // Not BOINC
	if ( main_exit_callback ) {
		main_exit_callback();
		std::exit( status );
	} else {
#ifndef _WIN32
		assert( false ); // Force a core dump for post-mortem debugging
#endif // Not _WIN32
		std::exit( status );
	}
#endif // BOINC

	throw UtilityExitException(file, line,  message);

}


/// @brief Conditional Exit with file + line + message + optional status
int
cond_exit(
	bool condition,
	char const * file,
	int const line,
	std::string const & message,
	int const status
){
	if ( condition ) return 1;
	exit( file, line, message, status );
	return 0; // keep compiler happy.
}


} // namespace utility
