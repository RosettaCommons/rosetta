// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/CarbohydrateInfoManager.cc
/// @brief   Method definitions for CarbohydrateInfoManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>

// Utility header
#include <utility/exit.hh>

// Basic header
#include <basic/database/open.hh>

// C++ header
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


// Singleton set-up
namespace utility {

using core::chemical::carbohydrates::CarbohydrateInfoManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< CarbohydrateInfoManager >::singleton_mutex_ {};
template <> std::atomic< CarbohydrateInfoManager * > utility::SingletonBase< CarbohydrateInfoManager >::instance_( 0 );
#else
template <> CarbohydrateInfoManager * utility::SingletonBase< CarbohydrateInfoManager >::instance_( 0 );
#endif

}  // namespace utility


namespace core {
namespace chemical {
namespace carbohydrates {

using namespace core;


// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
bool
CarbohydrateInfoManager::is_valid_sugar_code( std::string const & code )
{
	return get_instance()->code_to_root_map().count( code );
}

std::string const &
CarbohydrateInfoManager::root_from_code( std::string const & code )
{
	if ( ! is_valid_sugar_code( code ) ) {
		utility_exit_with_message( "CarbohydrateInfoManager::root_from_code(): " +
			code + " is not a valid 3-letter code for carbohydrates.  " +
			"Has it been added to database/chemical/carbohydrates/codes_to_roots.map?");
	}
	return get_instance()->code_to_root_map().find( code )->second;
}


bool
CarbohydrateInfoManager::is_valid_ring_affix( char affix )
{
	return get_instance()->ring_affixes().contains( affix );
}

char
CarbohydrateInfoManager::ring_affix_from_ring_size( core::Size ring_size )
{
	return get_instance()->ring_size_to_morphemes_map().find( ring_size )->second.first;
}

std::string const &
CarbohydrateInfoManager::morpheme_from_ring_size( core::Size ring_size )
{
	return get_instance()->ring_size_to_morphemes_map().find( ring_size )->second.second;
}


bool
CarbohydrateInfoManager::is_valid_modification_affix( std::string const & affix )
{
	return get_instance()->affix_to_patch_map().count( affix );
}

std::string const &
CarbohydrateInfoManager::patch_name_from_affix( std::string const & affix )
{
	return get_instance()->affix_to_patch_map().find( affix )->second;
}

/// @return  A position from 1 to 9 or 0 if there is no default value.
core::uint
CarbohydrateInfoManager::default_position_from_affix( std::string const & affix )
{
	return get_instance()->affix_to_position_map().find( affix )->second;
}


VariantType
CarbohydrateInfoManager::branch_variant_type_from_position( core::uint position )
{
	switch ( position ) {
		case 1:
			return C1_BRANCH_POINT;
		case 2:
			return C2_BRANCH_POINT;
		case 3:
			return C3_BRANCH_POINT;
		case 4:
			return C4_BRANCH_POINT;
		case 5:
			return C5_BRANCH_POINT;
		case 6:
			return C6_BRANCH_POINT;
		case 7:
			return C7_BRANCH_POINT;
		case 8:
			return C8_BRANCH_POINT;
		case 9:
			return C9_BRANCH_POINT;
		default:
			return NO_VARIANT;
	}
}


// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
CarbohydrateInfoManager::CarbohydrateInfoManager() {}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
CarbohydrateInfoManager *
CarbohydrateInfoManager::create_singleton_instance()
{
	return new CarbohydrateInfoManager;
}


// Get the map of Rosetta PDB 3-letter codes for saccharide residues mapped to the corresponding root requested,
// creating them if necessary.
// Called by the public static method root_from_code().
std::map< std::string, std::string > const &
CarbohydrateInfoManager::code_to_root_map()
{
	// Only create map one time, as needed.
	if ( code_to_root_map_.empty() ) {
		code_to_root_map_ = read_codes_and_roots_from_database_file(
				basic::database::full_name( "chemical/carbohydrates/codes_to_roots.map" ) );
	}
	return code_to_root_map_;
}


// Get the map of carbohydrate ring sizes and their 1-letter affixes and morphemes requested, creating it if
// necessary.
// Called by the public static methods ring_affix_from_ring_size() and morpheme_from_ring_size().
std::map< core::Size, std::pair< char, std::string > > const &
CarbohydrateInfoManager::ring_size_to_morphemes_map()
{
	// Only create map one time, as needed.
	if ( ring_size_to_morphemes_map_.empty() ) {
		ring_size_to_morphemes_map_ = read_ring_sizes_and_morphemes_from_database_file(
				basic::database::full_name( "chemical/carbohydrates/ring_size_to_morphemes.map" ) );
	}
	return ring_size_to_morphemes_map_;
}

// Get a list of valid 1-letter affixes for ring size, creating it if necessary.
utility::vector1< char > const &
CarbohydrateInfoManager::ring_affixes()
{
	using namespace std;
	typedef map< Size, pair< char, string > > Map;

	// Only create list one time, as needed.
	if ( ring_affixes_.empty() ) {
		Map const & map( ring_size_to_morphemes_map() );
		for ( Map::const_iterator it( map.begin() ), it_end( map.end() ); it != it_end; ++it ) {
			char affix( it->second.first );
			if ( affix != '\0' ) {
				ring_affixes_.push_back( affix );
			}
		}
	}
	return ring_affixes_;
}


// Get the table of nomenclature data for sugar modifications, creating it if necessary.
SugarModificationsNomenclatureTable const &
CarbohydrateInfoManager::nomenclature_table()
{
	// Only create table one time, as needed.
	if ( nomenclature_table_.empty() ) {
		nomenclature_table_ = read_nomenclature_table_from_database_file(
				basic::database::full_name( "chemical/carbohydrates/sugar_modifications.table" ) );
	}
	return nomenclature_table_;
}

// Get a map of sugar modification affixes to Rosetta patch names, creating it if necessary.
std::map< std::string, std::string > const &
CarbohydrateInfoManager::affix_to_patch_map()
{
	using namespace std;

	// Only create map one time, as needed.
	if ( affix_to_patch_map_.empty() ) {
		SugarModificationsNomenclatureTable const & table( nomenclature_table() );
		for ( SugarModificationsNomenclatureTable::const_iterator it( table.begin() ), it_end( table.end() );
				it != it_end; ++it ) {
			affix_to_patch_map_[ it->second.short_affix ] = it->second.patch_name;
		}
	}
	return affix_to_patch_map_;
}

// Get a map of sugar modification affixes to default positions, creating it if necessary.
std::map< std::string, core::uint > const &
CarbohydrateInfoManager::affix_to_position_map()
{
	using namespace std;

	// Only create map one time, as needed.
	if ( affix_to_position_map_.empty() ) {
		SugarModificationsNomenclatureTable const & table( nomenclature_table() );
		for ( SugarModificationsNomenclatureTable::const_iterator it( table.begin() ), it_end( table.end() );
				it != it_end; ++it ) {
			affix_to_position_map_[ it->second.short_affix ] = it->second.default_position;
		}
	}
	return affix_to_position_map_;
}

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
