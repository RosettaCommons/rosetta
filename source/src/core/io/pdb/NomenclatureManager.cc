// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/RingConformerManager.hh
/// @brief   Method definitions for NomenclatureManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/io/pdb/NomenclatureManager.hh>

// Package header
#include <core/io/pdb/alt_codes_io.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/database/open.hh>

// C++ headers
#include <string>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


// Singleton set-up
namespace utility {

using core::io::pdb::NomenclatureManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< NomenclatureManager >::singleton_mutex_ {};
template <> std::atomic< NomenclatureManager * > utility::SingletonBase< NomenclatureManager >::instance_( 0 );
#else
template <> NomenclatureManager * utility::SingletonBase< NomenclatureManager >::instance_( 0 );
#endif

}  // namespace utility


// Construct tracer.
static thread_local basic::Tracer TR( "core.io.pdb.NomenclatureManager" );


namespace core {
namespace io {
namespace pdb {

// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
/// @return  If a file of alternative 3-letter codes has been provided at the command line, return a pair of Rosetta
/// names (3-letter code and base ResidueType name, if available); if no file has been provided, or if there is no
/// alternative for the given code, simply return the PDB code and an empty string.
std::pair< std::string, std::string >
NomenclatureManager::rosetta_names_from_pdb_code( std::string const & pdb_code )
{
	using namespace std;
	using namespace basic::options;

	if ( option[ OptionKeys::in::alternate_3_letter_codes ].active() ) {  // Are alternate codes allowed?
		AltCodeMap const & alt_codes( get_instance()->get_alternate_3_letter_code_map() );
		AltCodeMap::const_iterator alt_code_pair( alt_codes.find( pdb_code ) );
		if ( alt_code_pair != alt_codes.end() ) {  // Is there an alternate for this code?
			// Get the value of this key/value pair.
			pair< string, string > const & rosetta_names( alt_code_pair->second );
			TR << "Accepting alternate code " << rosetta_names.first << " for " << pdb_code << '.' << endl;
			return rosetta_names;
		}
	}
	return std::make_pair( pdb_code, "" );
}


// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
NomenclatureManager::NomenclatureManager() {}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
NomenclatureManager *
NomenclatureManager::create_singleton_instance()
{
	return new NomenclatureManager;
}


// Get the map requested, creating it if necessary.
// Called by the public static method rosetta_names_from_pdb_code()
AltCodeMap const &
NomenclatureManager::get_alternate_3_letter_code_map()
{
	using namespace std;
	using namespace basic::options;
	using namespace core;

	// Only create map one time, as needed.
	if ( alt_codes_.size() == 0 ) {
		utility::vector1< string > const & file_list( option[ OptionKeys::in::alternate_3_letter_codes ]() );
		Size const n_files( file_list.size() );
		for ( uint i( 1 ); i <= n_files; ++i ) {
			string const & filename( find_alternate_codes_file( file_list[ i ] ) );
			AltCodeMap alt_codes( read_alternative_3_letter_codes_from_database_file( filename ) );
			alt_codes_.insert( alt_codes.begin(), alt_codes.end() );
		}
	}
	return alt_codes_;
}

// Try various combinations to locate the specific file being requested by the user.
// (inspired by core::scoring::ScoreFunction::find_weights_file())
std::string
NomenclatureManager::find_alternate_codes_file( std::string filename )
{
	using namespace utility::io;

	std::string const & path( basic::database::full_name( "input_output/3-letter_codes/" ) );
	std::string const ext( ".codes" );

	izstream potential_file( filename );
	if ( potential_file.good() ) {
		return filename;
	} else {
		izstream potential_file( filename + ext );  // Perhaps the user didn't use the .codes extension.
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
					utility_exit_with_message( "Unable to open alternative codes file. Neither ./" + filename +
							" nor " + "./" + filename + ext +
							" nor " + path + filename +
							" nor " + path + filename + ext + " exists." );
				}
			}
		}
	}
	return "INCONCEIVABLE!";  // Code can never reach here.
}

}  // namespace pdb
}  // namespace io
}  // namespace core
