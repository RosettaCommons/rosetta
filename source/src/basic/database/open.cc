// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/database/open.cc
/// @brief  Functions for opening database files
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Unit headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Project headers
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/PathName.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <cstdlib>
#include <iostream>
using basic::T;

//Auto Headers
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

namespace basic {
namespace database {

static THREAD_LOCAL basic::Tracer TR( "basic.io.database" );


/// @brief Open a database file on a provided stream
bool
open(
	utility::io::izstream & db_stream,
	std::string const & db_file,
	bool warn /* = true */
)
{
	using namespace utility::excn;

	if ( db_stream.good() ) {
		db_stream.close();
		db_stream.clear();
	}
	if ( db_file.length() == 0 ) {
		throw EXCN_Msg_Exception("Unable to open database file ''");
		return false;
	}

	db_stream.open( full_name( db_file, warn ) );

	if ( db_stream ) { // Open succeeded
		TR << "Database file opened: " << db_file << std::endl;
		return true;
	} else { // Open failed
		std::stringstream err_msg;
		err_msg
			<< "Database file open failed for: \"" << db_file << "\"" << std::endl;
		throw EXCN_Msg_Exception(err_msg.str());

#ifdef __native_client__
		throw( "ERROR: Database file open failed for: " + db_file );
#endif
		db_stream.close();
		db_stream.clear();
		return false;
	}
}


/// @brief Full-path database file name
std::string
full_name(
	std::string const & db_file,
	bool warn // = true
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	for ( size_t i = 1, i_end = option[ in::path::database ]().size(); i <= i_end; ++i ) {
		std::string fname = option[ in::path::database ](i).name() + db_file;
		if ( utility::file::file_exists(fname) || utility::file::file_exists(fname + ".gz") ) return fname;
	}
	// Don't exit -- sometimes caller wants to check if file exists (e.g. Dunbrack .bin file)
	//utility_exit_with_message("Unable to locate database file "+db_file);
	if ( warn ) Warning() << "Unable to locate database file " << db_file << std::endl;
	return option[ in::path::database ](1).name() + db_file;
}


} // namespace database
} // namespace basic
