// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/antibody/grafting/RegExManager.cc
/// @brief   Method implementations for RegExManager.
/// @author  Brian D. Weitzner <brian.weitzner@gmail.com>

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

// Unit headers
#include <protocols/antibody/grafting/regex_manager.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// C++ header

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


// Singleton set-up
namespace utility {

using protocols::antibody::grafting::RegExManager;

#ifdef MULTI_THREADED
template <> std::mutex utility::SingletonBase< RegExManager >::singleton_mutex_ {};
template <> std::atomic< RegExManager * > utility::SingletonBase< RegExManager >::instance_( 0 );
#else
template <> RegExManager * utility::SingletonBase< RegExManager >::instance_( 0 );
#endif

}  // namespace utility


namespace protocols {
namespace antibody {
namespace grafting {

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting");

// using namespace std;
// using namespace core;


// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
std::string RegExManager::H1_pattern() const { return H1_pattern_; }
std::string RegExManager::H3_pattern() const { return H3_pattern_; }

std::string RegExManager::L1_pattern() const { return L1_pattern_; }
std::string RegExManager::L3_pattern() const { return L3_pattern_; }

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
RegExManager::RegExManager() { load_regex_from_db(); }

// Singleton-creation function for use with utility::thread::threadsafe_singleton
RegExManager *
RegExManager::create_singleton_instance()
{
	return new RegExManager;
}

/// @brief Loads regex from database, easier for user modification this way
void RegExManager::load_regex_from_db() {

	// find the regex database file
	// expected format is:

	// # comment
	// H1:"REGEX_STRING"

	using std::string;

	string dir = "protocol_data/antibody/";
	string filename = "cdr_regex.txt";
	string path_to_file = basic::database::find_database_path( dir, filename );

	std::ifstream f(path_to_file);
	string line;

	TR << "Path to REGEX file is: " << path_to_file << std::endl;

	struct {
		string cdr_name;
		string & pattern;
	} R[] {
		{ "H1", H1_pattern_ },
		{ "H3", H3_pattern_ },
		{ "L1", L1_pattern_ },
		{ "L3", L3_pattern_ },
	};

	while ( std::getline( f, line ) ) {
		// load pattern from database file into data member
		for( auto & region : R ) {
			if( utility::startswith( line, region.cdr_name ) ) {
				region.pattern = utility::strip( utility::split( line )[ 2 ], '\"' );
				TR << "The detected " << region.cdr_name << " REGEX is: " << region.pattern << std::endl;
			}
		}
	}
}

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
