// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/exit.cc
/// @brief  Program exit functions and macros
/// @author David Kim (dekim@u.washington.edu)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)

#ifdef USEMPI
#include <mpi.h>
#endif
// Unit headers
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
// C++ headers
#include <cassert>
#include <cstdlib>
#include <iostream>

// Boinc headers
#ifdef BOINC
#include <utility/io/izstream.hh>
#include <boinc_api.h>
#include <string>
#endif


namespace utility {


static void (*exit_callback)(void) = 0;

void set_exit_callback( void (*my_callback)(void) )
{
	exit_callback = my_callback;
}


class EXCN_utility_exit : public excn::EXCN_Base { //noboday should be allowed to throw this... that's why its privately hidden in this modul...
public:
	EXCN_utility_exit( std::string const& msg, std::string const& file, int const line );
	virtual void show( std::ostream& ) const;
private:
	std::string const msg_;
	std::string const file_;
	int const line_;
};


EXCN_utility_exit::EXCN_utility_exit( std::string const& msg, std::string const& file, int const line ) :
	msg_( msg ),
	file_( file ),
	line_( line )
{}

void EXCN_utility_exit::show( std::ostream& os ) const {
	os << "\n\n[ERROR] EXCN_utility_exit has been thrown from: "
		 << file_ << " line: " << line_ << "\n";
	if ( ! msg_.empty() ) os << "ERROR: " << msg_ << "\n\n";
}


/// @brief Exit with file + line + message + optional status
void
exit(
	std::string const & file,
	int const line,
	std::string const & message,
	int const status
)
{
#ifdef EXIT_THROWS_EXCEPTION
	throw EXCN_utility_exit( message, file, line );
#endif
	if ( ! message.empty() ) std::cerr << std::endl << "ERROR: " << message << std::endl;
	std::cerr << "ERROR:: Exit from: " << file << " line: " << line << std::endl;
	std::cerr.flush();

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
	if( exit_callback ) exit_callback();
	else {
		#ifndef _WIN32
			assert( false ); // Force a core dump for post-mortem debugging
		#endif // Not _WIN32
		std::exit( status );
	}
#endif // BOINC
}


/// @brief Conditional Exit with file + line + message + optional status
int
cond_exit(
	bool condition,
	std::string const & file,
	int const line,
	std::string const & message,
	int const status
){
	if( condition ) return 1;
	exit( file, line, message, status );
	return 0; // keep compiler happy.
}



} // namespace utility
