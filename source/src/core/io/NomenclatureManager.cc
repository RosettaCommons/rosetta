// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/NomenclatureManager.hh
/// @brief   Method definitions for NomenclatureManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/io/NomenclatureManager.hh>

// Package header
#include <core/io/alt_codes_io.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/io/izstream.hh>
#include <utility/io/util.hh>
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

using core::io::NomenclatureManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< NomenclatureManager >::singleton_mutex_ {};
template <> std::atomic< NomenclatureManager * > utility::SingletonBase< NomenclatureManager >::instance_( 0 );
#else
template <> NomenclatureManager * utility::SingletonBase< NomenclatureManager >::instance_( 0 );
#endif

}  // namespace utility


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.io.pdb.NomenclatureManager" );


namespace core {
namespace io {

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
			TR << "Accepting alternate code " << pdb_code << " for " << rosetta_names.first << '.' << endl;
			return rosetta_names;
		}
	}
	return std::make_pair( pdb_code, "" );
}

bool NomenclatureManager::is_NA( std::string const & name3 ) {
	std::set< std::string > const & na_set( get_instance()->na_set() );
	return na_set.count( name3 );
}

bool NomenclatureManager::is_old_RNA( std::string const & name3 ) {
	std::set< std::string > const & old_rna_set( get_instance()->old_rna_set() );
	return old_rna_set.count( name3 );
}

bool NomenclatureManager::is_old_DNA( std::string const & name3 ) {
	std::set< std::string > const & old_dna_set( get_instance()->old_dna_set() );
	return old_dna_set.count( name3 );
}

bool NomenclatureManager::decide_is_d_aa( std::string const & name3 ) {
	std::set< std::string > const & d_aa_set( get_instance()->d_aa_set() );
	return d_aa_set.count( name3 );
}

bool NomenclatureManager::decide_is_l_aa( std::string const & name3 ) {
	std::set< std::string > const & l_aa_set( get_instance()->l_aa_set() );
	return l_aa_set.count( name3 );
}

bool NomenclatureManager::decide_is_known_achiral( std::string const & name3 ) {
	std::set< std::string > const & achiral_set( get_instance()->achiral_set() );
	return achiral_set.count( name3 );
}

bool NomenclatureManager::is_metal( std::string const & name3 ) {
	std::set< std::string > const & metal_set( get_instance()->metal_set() );
	return metal_set.count( name3 );
}

bool NomenclatureManager::is_sugar( std::string const & name3 ) {
	std::set< std::string > const & sugar_set( get_instance()->sugar_set() );
	return sugar_set.count( name3 );
}


// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
NomenclatureManager::NomenclatureManager()
{
	using namespace std;
	using namespace basic::options;
	using namespace core;

	if ( option[ OptionKeys::in::alternate_3_letter_codes ].user() ) {
		utility::vector1< string > const & file_list( option[ OptionKeys::in::alternate_3_letter_codes ]() );
		Size const n_files( file_list.size() );
		for ( uint i( 1 ); i <= n_files; ++i ) {
			string const & filename( find_alternate_codes_file( file_list[ i ] ) );
			AltCodeMap alt_codes( read_alternative_3_letter_codes_from_database_file( filename ) );
			alt_codes_.insert( alt_codes.begin(), alt_codes.end() );
		}
	}

	utility::vector1< string > const & file_list2 = option[ OptionKeys::in::name3_property_codes ]();
	for ( Size jj = 1; jj <= file_list2.size(); ++jj ) {
		utility::vector1< string > const lines( utility::io::get_lines_from_file_data( find_alternate_codes_file( file_list2[ jj ] ) ) );
		for ( Size ii = 1; ii <= lines.size(); ++ii ) {

			if ( lines[ ii ].size() == 0 || lines[ ii ][ 0 ] == '#' ) continue;
			istringstream word_by_word( lines[ ii ] );

			string name3, value;
			getline( word_by_word, value, '\t' );
			getline( word_by_word, name3, '\t' );

			if ( value == "NA" ) {
				is_NA_.insert( name3 );
			} else if ( value == "OLD_DNA" ) {
				is_old_DNA_.insert( name3 );
			} else if ( value == "OLD_RNA" ) {
				is_old_RNA_.insert( name3 );
			} else if ( value == "L_AA" ) {
				l_aa_set_.insert( name3 );
			} else if ( value == "D_AA" ) {
				d_aa_set_.insert( name3 );
			} else if ( value == "ACHIRAL" ) {
				achiral_set_.insert( name3 );
			} else if ( value == "METAL" ) {
				metal_set_.insert( name3 );
			} else if ( value == "SUGAR" ) {
				sugar_set_.insert( name3 );
			} else {
				utility_exit_with_message( "Line in name3 properties file \"" + file_list2[ jj ] + "\" malformed:\"" + lines[ ii ] + "\"" );
			}
		}
	}
}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
NomenclatureManager *
NomenclatureManager::create_singleton_instance()
{
	return new NomenclatureManager;
}


// Get the map requested, creating it if necessary.
// Called by the public static method rosetta_names_from_pdb_code()
AltCodeMap const &
NomenclatureManager::get_alternate_3_letter_code_map() const
{
	return alt_codes_;
}

// Try various combinations to locate the specific file being requested by the user.
// (inspired by core::scoring::ScoreFunction::find_weights_file())
std::string
NomenclatureManager::find_alternate_codes_file( std::string const & filename )
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


}  // namespace io
}  // namespace core
