// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/carbohydrates/CarbohydrateInfoManager.cc
/// @brief   Method definitions for CarbohydrateInfoManager.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Vikram K. Mulligan (vmullig@uw.edu) -- Made the CarbohydrateInfoManager threadsafe.

// Unit headers
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/io/util.hh>
#include <utility/exit.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <basic/database/open.hh>

// C++ header
#include <map>


namespace core {
namespace chemical {
namespace carbohydrates {

using namespace std;
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
	return get_instance()->code_to_root_map().find( code )->second.first;
}

/// @return  L, D, or *, where * indicates that the stereochemistry is inherent to the name.
char
CarbohydrateInfoManager::default_stereochem_from_code( std::string const & code )
{
	return get_instance()->code_to_root_map().find( code )->second.second;
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

bool
CarbohydrateInfoManager::affix_has_inherent_position( std::string const & affix )
{
	return get_instance()->affix_to_position_inherency_map().find( affix )->second;
}


bool
CarbohydrateInfoManager::pair_has_linkage_statistics( std::string const & res1, std::string const & res2 )
{
	pair< string, string > const key( convert_residue_names_into_linkage_map_key( res1, res2 ) );
	//std::cout << "Map Keys: " << key.first << " " << key.second << std::endl;
	return get_instance()->linkage_conformers_map().count( key );
}

utility::vector1< LinkageConformerData >
CarbohydrateInfoManager::linkages_from_pair( std::string const & res1, std::string const & res2 )
{
	pair< string, string > const key( convert_residue_names_into_linkage_map_key( res1, res2 ) );
	return get_instance()->linkage_conformers_map().find( key )->second;
}


VariantType
CarbohydrateInfoManager::branch_variant_type_from_atom_name( std::string const & atom_name )
{
	if ( atom_name == "O1" ) {
		return C1_BRANCH_POINT;
	} else if ( atom_name == "O2" ) {
		return C2_BRANCH_POINT;
	} else if ( atom_name == "O3" ) {
		return C3_BRANCH_POINT;
	} else if ( atom_name == "O4" ) {
		return C4_BRANCH_POINT;
	} else if ( atom_name == "O5" ) {
		return C5_BRANCH_POINT;
	} else if ( atom_name == "O6" ) {
		return C6_BRANCH_POINT;
	} else if ( atom_name == "O7" ) {
		return C7_BRANCH_POINT;
	} else if ( atom_name == "O8" ) {
		return C8_BRANCH_POINT;
	} else if ( atom_name == "O9" ) {
		return C9_BRANCH_POINT;
	} else {
		return NO_VARIANT;
	}
}

VariantType
CarbohydrateInfoManager::branch_variant_type_from_position( core::uint const position )
{
	switch ( position ) {
	case 1 :
		return C1_BRANCH_POINT;
	case 2 :
		return C2_BRANCH_POINT;
	case 3 :
		return C3_BRANCH_POINT;
	case 4 :
		return C4_BRANCH_POINT;
	case 5 :
		return C5_BRANCH_POINT;
	case 6 :
		return C6_BRANCH_POINT;
	case 7 :
		return C7_BRANCH_POINT;
	case 8 :
		return C8_BRANCH_POINT;
	case 9 :
		return C9_BRANCH_POINT;
	default :
		return NO_VARIANT;
	}
}


// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
CarbohydrateInfoManager::CarbohydrateInfoManager() {}

// Get the map of Rosetta PDB 3-letter codes for saccharide residues mapped to the corresponding root requested,
// creating them if necessary.
// Called by the public static method root_from_code().
std::map< std::string, std::pair< std::string, char > > const &
CarbohydrateInfoManager::code_to_root_map()
{
#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(code_to_root_map_mutex_);
		isempty = code_to_root_map_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(code_to_root_map_mutex_);
		if ( code_to_root_map_.empty() ) {
			code_to_root_map_ = read_codes_and_roots_from_database_file( basic::database::full_name( "chemical/carbohydrates/codes_to_roots.map" ) );
		}
	}
#else
	// Only create map one time, as needed.
	if ( code_to_root_map_.empty() ) {
		code_to_root_map_ = read_codes_and_roots_from_database_file(
			basic::database::full_name( "chemical/carbohydrates/codes_to_roots.map" ) );
	}
#endif
	return code_to_root_map_;
}


// Get the map of carbohydrate ring sizes and their 1-letter affixes and morphemes requested, creating it if
// necessary.
// Called by the public static methods ring_affix_from_ring_size() and morpheme_from_ring_size().
std::map< core::Size, std::pair< char, std::string > > const &
CarbohydrateInfoManager::ring_size_to_morphemes_map()
{

#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(ring_size_to_morphemes_mutex_);
		isempty = ring_size_to_morphemes_map_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(ring_size_to_morphemes_mutex_);
		if ( ring_size_to_morphemes_map_.empty() ) {
			ring_size_to_morphemes_map_ = read_ring_sizes_and_morphemes_from_database_file( basic::database::full_name( "chemical/carbohydrates/ring_size_to_morphemes.map" ) );
		}
	}
#else
	// Only create map one time, as needed.
	if ( ring_size_to_morphemes_map_.empty() ) {
		ring_size_to_morphemes_map_ = read_ring_sizes_and_morphemes_from_database_file(
			basic::database::full_name( "chemical/carbohydrates/ring_size_to_morphemes.map" ) );
	}
#endif
	return ring_size_to_morphemes_map_;
}

// Get a list of valid 1-letter affixes for ring size, creating it if necessary.
utility::vector1< char > const &
CarbohydrateInfoManager::ring_affixes()
{
	typedef map< Size, pair< char, string > > Map;
#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(ring_affixes_mutex_);
		isempty = ring_affixes_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(ring_affixes_mutex_);
		if ( ring_affixes_.empty() ) {
			Map const & map( ring_size_to_morphemes_map() );
			for ( Map::const_iterator it( map.begin() ), it_end( map.end() ); it != it_end; ++it ) {
				char affix( it->second.first );
				if ( affix != '\0' ) {
					ring_affixes_.push_back( affix );
				}
			}
		}
	}
#else
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
#endif
	return ring_affixes_;
}


// Get the table of nomenclature data for sugar modifications, creating it if necessary.
SugarModificationsNomenclatureTable const &
CarbohydrateInfoManager::nomenclature_table()
{
#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(nomenclature_table_mutex_);
		isempty = nomenclature_table_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(nomenclature_table_mutex_);
		if ( nomenclature_table_.empty() ) {
			nomenclature_table_ = read_nomenclature_table_from_database_file(
				basic::database::full_name( "chemical/carbohydrates/sugar_modifications.table" ) );
		}
	}
#else
	// Only create table one time, as needed.
	if ( nomenclature_table_.empty() ) {
		nomenclature_table_ = read_nomenclature_table_from_database_file(
			basic::database::full_name( "chemical/carbohydrates/sugar_modifications.table" ) );
	}
#endif
	return nomenclature_table_;
}

// Get a map of sugar modification affixes to Rosetta patch names, creating it if necessary.
std::map< std::string, std::string > const &
CarbohydrateInfoManager::affix_to_patch_map()
{
#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(affix_to_patch_mutex_);
		isempty = affix_to_patch_map_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(affix_to_patch_mutex_);
		if ( affix_to_patch_map_.empty() ) {
			SugarModificationsNomenclatureTable const & table( nomenclature_table() );
			for ( SugarModificationsNomenclatureTable::const_iterator it( table.begin() ), it_end( table.end() );
					it != it_end; ++it ) {
				affix_to_patch_map_[ it->second.short_affix ] = it->second.patch_name;
			}
		}
	}
#else
	// Only create map one time, as needed.
	if ( affix_to_patch_map_.empty() ) {
		SugarModificationsNomenclatureTable const & table( nomenclature_table() );
		for ( SugarModificationsNomenclatureTable::const_iterator it( table.begin() ), it_end( table.end() );
				it != it_end; ++it ) {
			affix_to_patch_map_[ it->second.short_affix ] = it->second.patch_name;
		}
	}
#endif
	return affix_to_patch_map_;
}

// Get a map of sugar modification affixes to default positions, creating it if necessary.
std::map< std::string, core::uint > const &
CarbohydrateInfoManager::affix_to_position_map()
{
#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(affix_to_position_mutex_);
		isempty = affix_to_position_map_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(affix_to_position_mutex_);
		// Only create map one time, as needed.
		if ( affix_to_position_map_.empty() ) {
			SugarModificationsNomenclatureTable const & table( nomenclature_table() );
			for ( SugarModificationsNomenclatureTable::const_iterator it( table.begin() ), it_end( table.end() );
					it != it_end; ++it ) {
				affix_to_position_map_[ it->second.short_affix ] = it->second.default_position;
			}
		}
	}
#else
	// Only create map one time, as needed.
	if ( affix_to_position_map_.empty() ) {
		SugarModificationsNomenclatureTable const & table( nomenclature_table() );
		for ( SugarModificationsNomenclatureTable::const_iterator it( table.begin() ), it_end( table.end() );
				it != it_end; ++it ) {
			affix_to_position_map_[ it->second.short_affix ] = it->second.default_position;
		}
	}
#endif
	return affix_to_position_map_;
}

// Get a map of sugar modification affixes to a boolean indication of if the position is inherent,
// creating it if necessary.
std::map< std::string, bool > const &
CarbohydrateInfoManager::affix_to_position_inherency_map()
{
#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(affix_to_position_inherency_mutex_);
		isempty = affix_to_position_inherency_map_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(affix_to_position_inherency_mutex_);
		// Only create map one time, as needed.
		if ( affix_to_position_inherency_map_.empty() ) {
			SugarModificationsNomenclatureTable const & table( nomenclature_table() );
			for ( SugarModificationsNomenclatureTable::const_iterator it( table.begin() ), it_end( table.end() );
					it != it_end; ++it ) {
				affix_to_position_inherency_map_[ it->second.short_affix ] = it->second.has_inherent_position;
			}
		}
	}
#else
	// Only create map one time, as needed.
	if ( affix_to_position_inherency_map_.empty() ) {
		SugarModificationsNomenclatureTable const & table( nomenclature_table() );
		for ( SugarModificationsNomenclatureTable::const_iterator it( table.begin() ), it_end( table.end() );
				it != it_end; ++it ) {
			affix_to_position_inherency_map_[ it->second.short_affix ] = it->second.has_inherent_position;
		}
	}
#endif
	return affix_to_position_inherency_map_;
}


// Get a map of linkage conformer statistical data, creating it if necessary.
LinkageConformers const &
CarbohydrateInfoManager::linkage_conformers_map()
{
	using namespace basic::options;
#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(linkage_conformers_mutex_);
		isempty = linkage_conformers_map_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(linkage_conformers_mutex_);
		// Only create map one time, as needed.
		if ( linkage_conformers_map_.empty() ) {
			string const filename( find_linkage_conformer_data_file(
				option[ OptionKeys::carbohydrates::linkage_conformer_data_file ]() ) );  // default is "default.table"
			linkage_conformers_map_ = read_linkage_conformers_from_database_file( filename );
		}
	}
#else
	// Only create map one time, as needed.
	if ( linkage_conformers_map_.empty() ) {
		string const filename( find_linkage_conformer_data_file(
			option[ OptionKeys::carbohydrates::linkage_conformer_data_file ]() ) );  // default is "default.table"
		linkage_conformers_map_ = read_linkage_conformers_from_database_file( filename );
	}
#endif
	return linkage_conformers_map_;
}

std::map< std::string, std::string > const &
CarbohydrateInfoManager::get_short_name_to_iupac_strings_map(){

	return get_instance()->short_name_to_iupac_strings_map();


}

std::map< std::string, std::string > const &
CarbohydrateInfoManager::short_name_to_iupac_strings_map(){
#ifdef MULTI_THREADED
	bool isempty(false);
	{
		utility::thread::ReadLockGuard readlock(short_name_to_iupac_strings_mutex_);
		isempty = short_name_to_iupac_strings_map_.empty();
	}

	if ( isempty ) {
		utility::thread::WriteLockGuard writelock(short_name_to_iupac_strings_mutex_);
		if ( short_name_to_iupac_strings_map_.empty() ) {
			std::string dir = "chemical/carbohydrates/common_glycans/";
			std::string filename = "common_names";
			std::string ext = ".txt";
			std::string path = basic::database::find_database_path(dir, filename, ext);
			short_name_to_iupac_strings_map_ = read_short_names_to_iupac_format_string( dir, path );

		}
	}
#else
	if ( short_name_to_iupac_strings_map_.empty() ) {
		std::string dir = "chemical/carbohydrates/common_glycans/";
		std::string filename = "common_names";
		std::string ext = ".txt";
		std::string path = basic::database::find_database_path(dir, filename, ext);
		short_name_to_iupac_strings_map_ = read_short_names_to_iupac_format_string( dir, path );

	}
#endif

	return short_name_to_iupac_strings_map_;
}


// Try various combinations to locate the specific file being requested by the user.
// (inspired by core::scoring::ScoreFunction::find_weights_file())
std::string
CarbohydrateInfoManager::find_linkage_conformer_data_file( std::string filename )
{
	using namespace utility::io;

	std::string dir = "chemical/carbohydrates/linkage_conformers/";
	std::string ext = ".table";

	return basic::database::find_database_path( dir, filename, ext);
}


// Helper function ////////////////////////////////////////////////////////////
std::pair< std::string, std::string >
convert_residue_names_into_linkage_map_key( std::string const & name1, std::string const & name2 )
{
	string fixed_name1( name1 );
	string fixed_name2( name2 );

	// First, we remove any trailing hyphens.
	fixed_name1 = fixed_name1.erase( fixed_name1.find_last_not_of( "-" ) + 1 );
	fixed_name2 = fixed_name2.erase( fixed_name2.find_last_not_of( "-" ) + 1 );

	// Next, we do not need/want the anomeric designation on the reducing-end residue.
	if ( fixed_name1.size() > 9 ) {  // Names shorter than this are non-saccharides.
		if ( fixed_name1.substr( 5, 5 ) == "alpha" ) {
			fixed_name1.erase( 5, 6 );  // We have to erase the hyphen after "alpha" too.
		}
		if ( fixed_name1.substr( 5, 4 ) == "beta" ) {
			fixed_name1.erase( 5, 5 );  // We have to erase the hyphen after "beta" too.
		}
	}

	// Next, we do not need the main-chain connectivity of the non-reducing-end residue.
	if ( fixed_name2.substr( 0, 2 ) == "->" ) {  // We assume that this is the start of a "->n)-".
		fixed_name2.erase( 0, 5 );
	}

	// Finally, make a pair and return.
	return make_pair( fixed_name2, fixed_name1 );
}

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
