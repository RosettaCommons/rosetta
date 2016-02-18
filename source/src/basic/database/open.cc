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

#if defined(MAC) || defined(__APPLE__)  ||  defined(__OSX__) || defined(linux) || defined(__linux__) || defined(__linux)
// POSIX specific headers
#include <pwd.h>
#endif

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

/// @brief Find a path to a file.
///
/// Try various combinations to locate the specific file being requested by the user.
/// (inspired by core::scoring::ScoreFunction::find_weights_file())
///
/// @athor Labonte <JWLabonte@jhu.edu>
std::string
find_database_path( std::string dir, std::string filename)
{
	using namespace utility::io;

	std::string const & path( basic::database::full_name( dir ) );

	izstream potential_file( filename );
	if ( potential_file.good() ) {
		return filename;
	} else {
		izstream potential_file( path + filename);  // Let's assume it's in the database in the usual spot.
		if ( potential_file.good() ) {
			return path + filename;
		} else {
			utility_exit_with_message( "Unable to open file. Neither ./" + filename +
				" nor " + "./" + filename +
				" nor " + path + filename + " exists." );
		}
	}
	return "WHAT THE @#$%!";  // Code can never reach here.
}


/// @brief Find a path to a file.
///
/// Try various combinations to locate the specific file being requested by the user.
/// (inspired by core::scoring::ScoreFunction::find_weights_file())
///
/// @athor Labonte <JWLabonte@jhu.edu>
std::string
find_database_path( std::string dir, std::string filename, std::string ext)
{
	using namespace utility::io;

	std::string const & path( basic::database::full_name( dir ) );

	izstream potential_file( filename );
	if ( potential_file.good() ) {
		return filename;
	} else {
		izstream potential_file( filename + ext );  // Perhaps the user didn't use the .table extension.
		if ( potential_file.good() ) {
			return filename + ext;
		} else {
			izstream potential_file( path + filename);  // Let's assume it's in the database in the usual spot.
			if ( potential_file.good() ) {
				return path + filename;
			} else {
				izstream potential_file( path + filename + ext );  // last try
				if ( potential_file.good() ) {
					return path + filename + ext;
				} else {
					utility_exit_with_message( "Unable to open file. Neither ./" + filename +
						" nor " + "./" + filename + ext +
						" nor " + path + filename +
						" nor " + path + filename + ext + " exists." );
				}
			}
		}
	}
	return "WHAT THE @#$%!";  // Code can never reach here.
}


/// @brief Does cache file (absolute path) exist?
/// if dir_only is true, will return true if the cache file could be created.
bool
find_cache_file(
	std::string const & cache_file,
	bool dir_only
) {
	if ( ! dir_only ) {
		bool exists = (utility::file::file_exists(cache_file) || utility::file::file_exists(cache_file + ".gz") );
		if ( TR.Debug.visible() && exists ) {
			TR.Debug << "Using '" << cache_file << "' as the cached file." << std::endl;
		}
		return exists;
	} else {
		// Does the directory exist/can it be created?
		std::string cache_dir = utility::file::FileName( cache_file ).path();
		if ( ! utility::file::create_directory_recursive( cache_dir ) ) {
			return false;
		}
		// Can we write a file in the directory?
		// Note that we *don't* want to try actually writing the actual file, due to race conditions.
		std::string tempfilename( utility::file::create_temp_filename( cache_dir, "writability_check" ) );
		std::ofstream tempfile( tempfilename.c_str() );
		bool usable = tempfile.good();
		tempfile.close();
		utility::file::file_delete(tempfilename); // Has internal file exist checks.
		if ( TR.Debug.visible() && usable ) {
			TR.Debug << "Using '" << cache_dir << "' as a cache directory." << std::endl;
		}
		return usable;
	}
}

/// @brief Get the (absolute) path to a given cached file.
/// If source_file is given, it's the full path to the source database file that's being cached.
/// If for_writing is true, will only check that the given file would be creatable.
/// Will return an empty string if it can't find a cache file.
std::string
full_cache_name(
	std::string const & short_name,
	std::string const & source_file,
	bool for_writing
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string cache_name;

	// First try the specified cache directories, if possible.
	if ( option[ in::path::database_cache_dir ].user() ) {
		cache_name = std::string(option[ in::path::database_cache_dir ]()) + short_name;
		if ( find_cache_file( cache_name, for_writing ) ) {
			return cache_name;
		}
	}

	char const * path = getenv("ROSETTA3_DBCACHE");
	if ( path && strlen(path) > 0 ) {
		cache_name = std::string(path) + "/" + short_name;
		if ( find_cache_file( cache_name, for_writing ) ) {
			return cache_name;
		}
	}

	// Then try the database directory
	// We don't iterate through all database directories, because in a multiple directory situation we don't want
	// to put the cache for one database into a different one
	if ( source_file.size() != 0 ) {
		cache_name = utility::file::FileName( source_file ).path() + utility::file::FileName( short_name ).bare_name();
		if ( find_cache_file( cache_name, for_writing ) ) {
			return cache_name;
		}
	}

	// No luck? Then fall back to the user's home directories.
	char const * homedir = getenv("XDG_CONFIG_HOME");
	if ( ! homedir || strlen(homedir) == 0 ) {
		homedir = getenv("HOME");
	}
#if defined(MAC) || defined(__APPLE__)  ||  defined(__OSX__) || defined(linux) || defined(__linux__) || defined(__linux)
	if ( ! homedir || strlen(homedir) == 0 ) {
		passwd const * unix_pwd( getpwuid(getuid()) );
		if ( unix_pwd ) {
			homedir = unix_pwd->pw_dir;
		}
	}
#endif

	if ( homedir && strlen(homedir) > 0 ) {
		cache_name = std::string(homedir) + "/.rosetta/database/" + short_name;
		if ( find_cache_file( cache_name, for_writing ) ) {
			return cache_name;
		}
	}

	return "";
}

} // namespace database
} // namespace basic
