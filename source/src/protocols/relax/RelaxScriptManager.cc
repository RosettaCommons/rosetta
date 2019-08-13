// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/relax/RelaxScriptManager.cc
/// @brief A singleton class for managing relax scripts, to ensure that they are loaded once and only once from disk.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <protocols/relax/RelaxScriptManager.hh>


// Unit headers

// Project header
#include <core/types.hh>

// Utility headers
#include <utility/pointer/memory.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/GeneralFileManager.hh>

// Core headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C++ headers
#include <string>
#include <fstream>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Construct tracer.
static basic::Tracer TR( "protocols.relax.RelaxScriptManager" );

namespace protocols {
namespace relax {


// EnergyMapContainer Public methods /////////////////////////////////////////////////////////////

/// @brief Constructor -- READS FROM DISK!
EnergyMapContainer::EnergyMapContainer( std::string const & filename ) :
	energy_map_( core::scoring::ScoreFunction::extract_weights_from_file( filename ) )
{}




// RelaxScriptFileContents Public methods /////////////////////////////////////////////////////////////
/// @brief File contents constructor.
RelaxScriptFileContents::RelaxScriptFileContents( utility::vector1< std::string > const & file_lines_in ) :
	utility::pointer::ReferenceCount(),
	file_lines_( file_lines_in )
{}

/// @brief Destructor.
RelaxScriptFileContents::~RelaxScriptFileContents() {}

/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
RelaxScriptFileContentsOP
RelaxScriptFileContents::clone() const {
	return utility::pointer::make_shared< RelaxScriptFileContents >(*this);
}




// RelaxScriptManager Public methods /////////////////////////////////////////////////////////////
// Static constant data access

/// @brief Get a relax script.  Load it from disk if it has not already been loaded.
/// @details Threadsafe and lazily loaded.  Requires the FastRelax or FastDesign mover's scorefunction to be provided.
RelaxScriptFileContents const &
RelaxScriptManager::get_relax_script(
	std::string const & filename,
	core::scoring::ScoreFunctionCOP const & mover_sfxn,
	bool const dualspace
) const {
	initialize_relax_scripts_in_database(); //Does nothing if already initialized.

	core::scoring::ScoreFunctionCOP sfxn( mover_sfxn == nullptr ? core::scoring::get_score_function( true ) : mover_sfxn );
	std::string const nearest_sfxn( get_nearest_sfxn_if_in_database( filename, sfxn ) );

	boost::function< RelaxScriptFileContentsOP () > creator( boost::bind( &RelaxScriptManager::create_relax_script_instance, this, boost::cref( filename ), nearest_sfxn, dualspace ) );
	return *( utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( relax_script_mutex_ ), std::make_tuple( filename, nearest_sfxn, dualspace ), filename_to_filecontents_map_ ) );
}

// RelaxScriptManager Private methods ////////////////////////////////////////////////////////////

/// @brief Empty constructor.
RelaxScriptManager::RelaxScriptManager() :
	SingletonBase< RelaxScriptManager >(),
	filename_to_filecontents_map_()
#ifdef MULTI_THREADED
	,
	relax_script_mutex_()
#endif //MULTI_THREADED
{}

/// @brief Create an instance of a RelaxScriptFileContents object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of RelaxScriptManager.
RelaxScriptFileContentsOP
RelaxScriptManager::create_relax_script_instance(
	std::string const & filename,
	std::string const & nearest_sfxn_name,
	bool const dualspace
) const {
	utility::vector1< std::string > filelines;

	std::string const filename_full( script_file_path( filename, nearest_sfxn_name, dualspace ) );

	std::ifstream infile( filename_full.c_str() );
	if ( !infile.good() ) {
		utility_exit_with_message( "[ERROR] Error opening relaxscript file '" + filename_full );
	}
	TR << "================== Reading script file: " << filename_full << " ==================" << std::endl;

	std::string line;
	while ( getline( infile, line ) ) {
		//store line
		filelines.push_back( line );
	}
	infile.close();

	for ( auto const & fl : filelines ) {
		TR << fl << std::endl;
	}

	return utility::pointer::make_shared< RelaxScriptFileContents >(filelines);
}


/// @brief Check whether (a) there is a period in the script file name, and, if so, (b) whether the base name exists in the
/// database.  If (a) is true, throw an error message, where the text depends on whether (b) is true.
/// @details PERFORMS DISK READ.
void
RelaxScriptManager::handle_script_file_with_period(
	std::string const & script_file,
	std::string::size_type position_of_first_period
){
	//In order to make a more helpful error message, check to see if basename exists in the database
	std::string const basename = script_file.substr( 0, position_of_first_period );
	std::string const basename_with_default_extension = basename + ".txt";

	bool const file_exists_with_default_extension =
		utility::file::file_exists( basic::database::full_name( "sampling/relax_scripts/" ) + basename_with_default_extension );
	if ( file_exists_with_default_extension ) {
		utility_exit_with_message(
			"[ERROR] relaxscript argument should not have extensions. "
			"Did you mean " + basename + " instead of " + script_file + "?" );
	} else {
		utility_exit_with_message(
			"[ERROR] relaxscript argument " + script_file +
			" should not have extensions. Additionally, " +
			basename + " does not appear to be a valid script name. "
			"Please look at main/database/sampling/relax_scripts/ or the wiki for valid names." );
	}
}

/// @brief Get the full script file path from the script name.
/// @details PEFRORMS DISK READ.
std::string
RelaxScriptManager::script_file_path(
	std::string const & script_file,
	std::string const & nearest_sfxn_name,
	bool const dualspace
) const {

	//General Case: First check if file exists locally, then check database. Add dualspace extension if necessary
	if ( utility::file::file_exists( script_file ) ) {
		//the user is providing their own script
		return script_file;
	} else { //let's look in the database
		//1. Make sure the name has no extension (because we will provide the correct extension here)
		std::string::size_type const position_of_first_period = script_file.find( '.' );
		if ( position_of_first_period != std::string::npos ) {
			handle_script_file_with_period( script_file, position_of_first_period );
		}

		std::string script_file_path(script_file);
		if ( nearest_sfxn_name != "" ) {
			script_file_path += "." + nearest_sfxn_name;
		}
		if ( dualspace ) {
			script_file_path += ".dualspace";
		}
		script_file_path += ".txt";

		TR << "Looking for " << script_file_path << std::endl;

		script_file_path = basic::database::find_database_path( "sampling/relax_scripts/", script_file_path, true );

		if ( script_file_path.empty() ) {
			std::string error_message = "Error! Could not find the following relax scripts in the database: " + script_file_path;
			utility_exit_with_message( error_message );
		}
		return script_file_path;
	}
	return script_file; //To keep compilers happy.
}

/// @brief Get the database scorefunction that is closest to the current one, if the relax script is a
/// database relax script that specifies a scorefunction.
std::string
RelaxScriptManager::get_nearest_sfxn_if_in_database(
	std::string const & script_file,
	core::scoring::ScoreFunctionCOP mover_sfxn
) const {
#ifdef MULTI_THREADED
	//Read lock for reading from relax_scripts_in_database_ (released at end of scope):
	utility::thread::ReadLockGuard readlock( relax_scripts_in_database_mutex_ );
#endif

	std::string returnstring;

	for ( DatabaseRelaxScript const & database_script : relax_scripts_in_database_ ) {
		if ( database_script.title == script_file ) {
			if ( database_script.sfxn_specialized ) {
				core::scoring::ScoreFunctionCOP sfxn( mover_sfxn );
				if ( sfxn == nullptr ) {
					sfxn = core::scoring::get_score_function( true );
				}
				returnstring = determine_closest_scorefunction( *sfxn, database_script.sfxns_supported );
			}
			break;
		}
	}
	return returnstring;
}

/// @brief Read the list of relax scripts from the database.
/// @details PEFRORMS DISK READ.
void
RelaxScriptManager::initialize_relax_scripts_in_database() const
{
#ifdef MULTI_THREADED
	//Check first with a read lock that we haven't already initialized:
	{ //Read-lock scope
		utility::thread::ReadLockGuard readlock( relax_scripts_in_database_mutex_ );
		if ( !relax_scripts_in_database_.empty() ) return;
	} //Release read lock and acquire write-lock:
	utility::thread::WriteLockGuard writelock( relax_scripts_in_database_mutex_ );
	//Check again that we haven't already initialized:
#endif
	if ( !relax_scripts_in_database_.empty() ) return;

	TR << "Reading relax scripts list from database." << std::endl;
	std::string const filename ( basic::database::find_database_path( "sampling/relax_scripts/", "index.dat", false ) );
	std::string const & index_file_contents ( utility::io::GeneralFileManager::get_instance()->get_file_contents( filename ) );
	std::vector< std::string > const index_file_lines( utility::split_by_newlines( index_file_contents ) );

	bool reached_title_row( false );
	std::string title, is_sfxn_dep, sfxns;

	for ( std::string const & line : index_file_lines ) {
		std::istringstream token_stream( line );
		token_stream >> title;
		if ( token_stream.fail() || title[0] == '#' ) continue;

		if ( ! reached_title_row ) {
			reached_title_row = true;
			continue;
		}

		token_stream >> is_sfxn_dep;
		if ( is_sfxn_dep == "Y" ) {
			token_stream >> sfxns;
			relax_scripts_in_database_.emplace_back( title, true, sfxns );
		} else {
			relax_scripts_in_database_.emplace_back( title, false );
		}
	}
}

/// @brief Get the difference between two EnergyMaps.
core::Real
RelaxScriptManager::sfxn_distance_squared(
	core::scoring::EnergyMap const & target_weights,
	core::scoring::EnergyMap const & candidate_weights
){
	core::Real sfxn_distance_squared = 0;
	for ( core::Size ii = core::scoring::fa_atr; ii <= core::scoring::n_score_types; ++ii ) {
		core::scoring::ScoreType ii_type = static_cast< core::scoring::ScoreType >( ii );
		if ( target_weights[ ii_type ] == 0.0 && candidate_weights[ ii_type ] == 0.0 ) continue;

		core::Real const difference = target_weights[ ii_type ] - candidate_weights[ ii_type ];
		sfxn_distance_squared += difference * difference;
	}

	return sfxn_distance_squared;
}

/// @brief Create an instance of an EnergyMapContainer.  Needed for threadsafe lazy loading.
/// @details Peforms onetime read from disk.
EnergyMapContainerCOP
RelaxScriptManager::create_energy_map_instance(
	std::string const & filename
) {
	std::string const full_filename( basic::database::full_name( "scoring/weights/" + filename + ".wts", true ) );
	return utility::pointer::make_shared< EnergyMapContainer >( full_filename );
}

/// @brief Given a scorefunction and a list of candidate scorefunctions, determine the closest one.
/// @details Performs onetime read from disk.
std::string
RelaxScriptManager::determine_closest_scorefunction(
	core::scoring::ScoreFunction const & target,
	utility::vector1< std::string > const & names_of_candidates
) const {
	core::scoring::EnergyMap const & target_weights = target.weights();

	bool first(true);
	core::Real best_sfxn_distance_squared( 0 );
	std::string best_candidate = "ref2015";

	for ( std::string const & candidate : names_of_candidates ) {

		/// @brief The following ensures that the energy maps are only loaded once, and cached:
		boost::function< EnergyMapContainerCOP () > creator( boost::bind( &RelaxScriptManager::create_energy_map_instance, boost::cref( candidate ) ) );
		core::scoring::EnergyMap const & dummy_weights( (utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( energy_maps_mutex_ ), candidate, energy_maps_map_ ) )->get_energy_map() );

		core::Real const score( sfxn_distance_squared( target_weights, dummy_weights ) );

		if ( first || score < best_sfxn_distance_squared ) {
			first = false;
			best_sfxn_distance_squared = score;
			best_candidate = candidate;
		}

		TR << candidate << " : " << score << std::endl;
	}

	return best_candidate;
}


} //protocols
} //relax
